library(SAJR)
library(biomaRt)
library(openxlsx)
library(dplyr)
ensembl=useMart("ensembl")
ensembl = useDataset("celegans_gene_ensembl",mart=ensembl)
gene_map=getBM(attributes=c("ensembl_gene_id", 'external_gene_name'), mart = ensembl, useCache = FALSE)
run_sajr <- function(grp1, grp1_samples, grp2, grp2_samples, gene_map){
  
  significance=0.05
  ann.gff='/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/sajr_output/at_pipe_test_asplicing.novel.gff'
  
  samples=c(grp1_samples,grp2_samples)
  cond=factor(c(rep(grp1,length(grp1_samples)),rep(grp2,length(grp2_samples))))
  cond=relevel(cond, ref=grp1)
  meta = list(cond=cond)
  meta
  
  # load annotation and count data
  data = loadSAData(ann.gff=ann.gff,samples)
  data<-setSplSiteTypes(data, ann.gff)
  
  # filter low values
  # new recommended filtering:
  data.f = data[data[['seg']][,'type'] %in% c('ALT','EXN','INT') & data[['seg']][,'position'] %in% c('LAST','INTERNAL','FIRST') & apply(data[['i']]+data[['e']]>=10,1,sum)>=2 & apply(data[['ir']],1,sd,na.rm=TRUE) > 0,]
  
  # test for significant differences
  data.f.glm = fitSAGLM(data.f,formula=terms(x ~ cond, keep.order=T),meta,0.05)
  data.f.pv = calcSAPvalue(data.f.glm)
  data.f[['seg']][,c("pvalue")]=data.f.pv[,"cond"]
  data.f[['seg']][,c("padj")] = p.adjust(data.f[['seg']][,c("pvalue")],method='BH')
  
  # calculate group means and log2FC
  data.f[['seg']][,c("psi.grp1")]<-rowMeans(data.f[['ir']][,1:length(grp1_samples), drop=FALSE], na.rm=T)
  data.f[['seg']][,c("psi.grp2")]<-rowMeans(data.f[['ir']][,(1+length(grp1_samples)):length(c(grp1_samples,grp2_samples)), drop=FALSE], na.rm=T)
  data.f[['seg']][,c("grp1.sd")]<-apply(data.f[['ir']][,1:length(grp1_samples), drop=FALSE],1,function(x) sd(x, na.rm=T))
  data.f[['seg']][,c("grp2.sd")]<-apply(data.f[['ir']][,(1+length(grp1_samples)):length(c(grp1_samples,grp2_samples)), drop=FALSE],1,function(x) sd(x, na.rm=T))
  data.f[['seg']][,c("psi.diff")]<-data.f[['seg']][,c("psi.grp2")] - data.f[['seg']][,c("psi.grp1")]
  data.f[['seg']][,c("log2FC")]<-log2((data.f[['seg']][,c("psi.grp2")]+0.01)/(data.f[['seg']][,c("psi.grp1")]+0.01)) # add one percent pseudoinclusion to avoid division-by-zero
  
  # read annotation data
  novelty<-read.delim("/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/sajr_output/novel_overlap_known_stringent_novel.at_pipe_test_asplicing.tsv", header=F)
  novelty<-novelty[!duplicated(novelty[,1]),]
  colnames(novelty)=c("segment","Novelty")
  rownames(novelty)=novelty[,'segment']
  
  names_from_overlap<-read.delim("/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/sajr_output/novel2known_from_overlap.at_pipe_test_asplicing.tsv", header=F)
  rownames(names_from_overlap)=names_from_overlap[,1]
  colnames(names_from_overlap)=c("sajrSegID","GeneName")
  head(names_from_overlap)
  
  names_from_comp<-read.delim("/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/sajr_output/sajr.at_pipe_test_asplicing.novel2known.tsv", header=F)
  rownames(names_from_comp)=names_from_comp[,1]
  colnames(names_from_comp) = c("gene_id","GeneNames","Class")
  head(names_from_comp)
  
  # save results
  data.res <- data.f
  colnames(data.res[['ir']]) <- paste("ir", colnames(data.res[['ir']]), sep = ".")
  colnames(data.res[['i']]) <- paste("i", colnames(data.res[['i']]), sep = ".")
  colnames(data.res[['e']]) <- paste("e", colnames(data.res[['e']]), sep = ".")
  res<-data.res[['seg']]
  res<-merge(res, data.res[['ir']], by="row.names")
  rownames(res)=res[,'Row.names']
  res<-res[,-1]
  res<-merge(res, data.res[['i']], by="row.names")
  rownames(res)=res[,'Row.names']
  res<-res[,-1]
  res<-merge(res, data.res[['e']], by="row.names")
  rownames(res)=res[,'Row.names']
  res<-res[,-1]
  res<-merge(res, novelty, by="row.names")
  rownames(res)=res[,'Row.names']
  res<-res[,-1]
  res<-merge(names_from_comp,res, by="gene_id")
  
  # sort the data according to p-value
  res<-res[order(res[,'pvalue']),]
  head(res)
  
  res[,'GeneNames'] <- gsub('gene:', '', res[,'GeneNames'])
  
  res = merge(res, gene_map, by.x = 'GeneNames', by.y = 'ensembl_gene_id', all.x = TRUE)
  # annotate genes by gene name
  
  write.table(res,file=paste(grp1, '_vs_', grp2, ".results.tsv", sep = ''), sep="\t", quote=F, row.names=F)
  write.xlsx(res, file = paste(grp1, '_vs_', grp2, ".results.xlsx", sep = ''))
  
}
# run SAJR pipeline for RNP6 samples

setwd('/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/sajr_output/count_files/')
samples = read.xlsx("/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/nextflow-alternative-splicing/sample_sheet.xlsx")
groups = unique(samples[,2])
for(g1 in 1:length(groups)){
  for(g2 in 1:length(groups)){
    if (g1 < g2) {
      print(paste('testing', groups[g1], 'and', groups[g2]))
      run_sajr(groups[g1], samples[samples[,2] == groups[g1],'sample_name'],  groups[g2],  samples[samples[,2] == groups[g2],'sample_name'], gene_map)
    }
  }
}
sink('/nexus/posix0/MAGE-flaski/service/hpc/home/sjiang/nextflow_asplicing_test/sajr_output/count_files/sessionInfo.txt' )
sessionInfo()
sink()