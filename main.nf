#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
      stageInMode 'symlink'
      stageOutMode 'move'

      script:
        """

        if [[ "${params.containers}" == "singularity" ]] ; 

          then

            cd ${params.image_folder}
        
            if [[ ! -f samtools-1.16.1.sif ]] ;
              then
                singularity pull samtools-1.16.1.sif docker://index.docker.io/mpgagebioinformatics/samtools:1.16.1
            fi
            
            if [[ ! -f rnaseq.python-3.8-1.sif ]] ;
              then
                singularity pull rnaseq.python-3.8-1.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-1
            fi

            if [[ ! -f sajr-1.0.0.sif ]] ;
              then
                singularity pull sajr-1.0.0.sif docker://index.docker.io/mpgagebioinformatics/sajr:1.0.0
            fi


        fi


        if [[ "${params.containers}" == "docker" ]] ; 

          then

            docker pull mpgagebioinformatics/samtools:1.16.1
            docker pull mpgagebioinformatics/rnaseq.python:3.8-1
            docker pull mpgagebioinformatics/sajr:1.0.0
        fi

        """
    }


process fill_config_template {
    stageInMode 'symlink'
    stageOutMode 'move'

    script: 

      """
      #!/usr/local/bin/python3
      
      # read in sample sheet
      # read in templates
      # replace place holders with file names
      
      import os
      import pandas as pd
      import openpyxl

      if os.path.exists("${params.tmp}"):
          pass
      else:
          os.makedirs("${params.tmp}")

      samplesheet = pd.read_excel("/workdir/nf-star-test/sample_sheet.xlsx", engine="openpyxl")
      samplesheet
      gff32sajr_template = open('${params.sajr_code}/sajr.gff32sajr.template.config')
      gff32sajr_config = gff32sajr_template.read()
      gff32sajr_template.close()
      # for sajr
      place_holder_bam_file='${params.sajr_output}${params.series}.merged.sorted.bam'
      place_holder_novel_gff='${params.sajr_output}${params.series}.novel.gff'
      place_holder_outfolder='${params.sajr_output}${params.series}/'
      place_holder_comp_out='${params.sajr_output}${params.series}.sajr.comp'
      # for gff32sajr
      place_holder_original_gff3='${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.nospace.gff3'
      place_holder_converted_gff='${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gff'
      place_holder_fasta_file='${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.fa'
      gff32sajr_config = gff32sajr_config.replace('[place_holder_original_gff3]', place_holder_original_gff3)
      gff32sajr_config = gff32sajr_config.replace('[place_holder_converted_gff]', place_holder_converted_gff)
      gff32sajr_config = gff32sajr_config.replace('[place_holder_fasta_file]', place_holder_fasta_file)
      gff32sajr_config = gff32sajr_config.replace('[place_holder_bam_file]', place_holder_bam_file)
      gff32sajr_config = gff32sajr_config.replace('[place_holder_novel_gff]', place_holder_novel_gff)
      gff32sajr_config = gff32sajr_config.replace('[place_holder_outfolder]', place_holder_outfolder)
      gff32sajr_config = gff32sajr_config.replace('[place_holder_comp_out]', place_holder_comp_out)
      gff32sajr_template = open('${params.tmp}/sajr.gff32sajr.${params.series}.config', 'w')
      gff32sajr_template.write(gff32sajr_config)
      gff32sajr_template.close()
      sajr_template = open('${params.sajr_code}/sajr.template.config')
      sajr_config = sajr_template.read()
      sajr_template.close()
      # iterate over samples and concatenate
      groups = dict()
      batch_in = list()
      batch_out = list()
      for index, row in samplesheet.iterrows():
          if not row[1] in groups:
              groups[row[1]] = 1
          else: 
              groups[row[1]] += 1     
          batch_in += [('${params.star_output}%s_%s.Aligned.sortedByCoord.out.bam' %(row[1], groups[row[1]]))]
          batch_out += [('${params.sajr_output}count_files/%s_%s' %(row[1], groups[row[1]]))]
      place_holder_batch_in=','.join(batch_in)
      place_holder_batch_out=','.join(batch_out)
      sajr_config = sajr_config.replace('[place_holder_batch_in]', place_holder_batch_in)
      sajr_config = sajr_config.replace('[place_holder_batch_out]', place_holder_batch_out)
      sajr_config = sajr_config.replace('[place_holder_bam_file]', place_holder_bam_file)
      sajr_config = sajr_config.replace('[place_holder_novel_gff]', place_holder_novel_gff)
      sajr_config = sajr_config.replace('[place_holder_outfolder]', place_holder_outfolder)
      sajr_config = sajr_config.replace('[place_holder_converted_gff]', place_holder_converted_gff)
      sajr_config = sajr_config.replace('[place_holder_fasta_file]', place_holder_fasta_file)
      sajr_config = sajr_config.replace('[place_holder_comp_out]', place_holder_comp_out)
      sajr_template = open('${params.tmp}/sajr.${params.series}.config', 'w')
      sajr_template.write(sajr_config)
      sajr_template.close()
      """
    
}

process sajr {
    stageInMode 'symlink'
    stageOutMode 'move'

    script: 

    """

    cd ${params.sajr_output}

    if [ ! -e ${params.organism}.${params.release}.gff ] ; then 

    # step 4
    # Convert a GTF reference to an SAJR specific GFF reference using SAJR’s annotation conversion mode.
    # sajr had an issue with spaces in the gff3 format
    
    sed 's/ /./g' ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gff3 > ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.nospace.gff3
   
    cd ${params.sajr_code}
    echo "step 4: convert gff3 to sajr gff"
    echo "java -jar sajr.jar gff32sajr ${params.tmp}/sajr.gff32sajr.${params.series}.config"
    java -jar sajr.jar gff32sajr ${params.tmp}/sajr.gff32sajr.${params.series}.config
    date
    
    fi

    ########################################

    if [ ! -s ${params.sajr_output}${params.series}.novel.gff ] ; then 
    # step 5
    # Run SAJR in de novo annotation mode to find novel splice-forms using the
    # merged BAM file and the known annotation to produce a novel annotation, novel.gff

    echo "step 5: find novel splice forms"
    echo "java -jar sajr.jar annotate ${params.tmp}/sajr.${params.series}.config" 
    cd ${params.sajr_code}
    java -jar sajr.jar annotate ${params.tmp}/sajr.${params.series}.config
    date
    fi

    ########################################

    if [ ! -s ${params.sajr_output}sajr.${params.series}.novel2known.tsv ] ; then 
    # step 6
    # Run SAJR in annotation comparison mode to compare the novel annotation with the known annotation and use
    # get_genename_from_junction_comparison.pl to filter the results:
    # get_genename_from_junction_comparison.pl sajr.comp > sajr.novel2known.tsv

    echo "step 6: compare novel splicing to known annotation"
    echo "java -jar sajr.jar sajrcomp ${params.tmp}/sajr.${params.series}.config"

    cd ${params.sajr_code}
    java -jar sajr.jar sajrcomp ${params.tmp}/sajr.${params.series}.config
    echo "perl get_genename_from_junction_comparison.pl ${params.sajr_output}${params.series}.sajr.comp > ${params.sajr_output}sajr.${params.series}.novel2known.tsv"
    perl get_genename_from_junction_comparison.pl ${params.sajr_output}${params.series}.sajr.comp > ${params.sajr_output}sajr.${params.series}.novel2known.tsv
    date

    fi

    ######################################## 

    if [ ! -s ${params.sajr_output}novel_overlap_known.${params.series}.gff ] ; then 
    # step 7
    # Use bedtools and get_genename_from_segment_overlap.pl to associate SAJR ids with known gene ids from the
    # reference: bedtools intersect -s -f 1.0 -loj -a novel.gff -b known.gff > novel_overlap_known.gff
   
    cd ${params.sajr_output}

    cp ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gff ./${params.organism}.${params.release}.gff
    echo "step 7: bedtools intersect novel and known isoforms"
    echo "bedtools intersect -s -f 1.0 -loj -a ${params.series}.novel.gff -b ${params.organism}.${params.release}.gff > novel_overlap_known.${params.series}.gff"
    bedtools intersect -s -f 1.0 -loj -a ${params.series}.novel.gff -b ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gff > novel_overlap_known.${params.series}.gff
    date
    fi

    ########################################

    if [ ! -s ${params.sajr_output}novel2known_from_overlap.${params.series}.tsv ] ; then 
    # step 8
    # get_genename_from_segment_overlap.pl novel_overlap_known.gff > novel2known_from_overlap.tsv
   
    cd ${params.sajr_code}
    echo "step 8: extract and add gene names"
    echo "perl get_genename_from_segment_overlap.pl ${params.sajr_output}novel_overlap_known.${params.series}.gff  > ${params.sajr_output}novel2known_from_overlap.${params.series}.tsv"
    perl get_genename_from_segment_overlap.pl ${params.sajr_output}novel_overlap_known.${params.series}.gff  > ${params.sajr_output}novel2known_from_overlap.${params.series}.tsv
    date
    fi

    ########################################

    if [ ! -s ${params.sajr_output}novel_overlap_known_stringent.${params.series}.gff ] ; then 
    # step 9
    # Use bedtools and annotate_novel_segments.pl to annotate novel spliced regions:
    # bedtools intersect -s -f 1.0 -r -loj -a novel.gff -b known.gff > novel_overlap_known_stringent.gff
   
    cd ${params.sajr_output}
    echo "step 9: bedtools stringent overlap with known genes"
    echo "bedtools intersect -s -f 1.0 -r -loj -a ${params.series}.novel.gff -b  ${params.organism}.${params.release}.gff > novel_overlap_known_stringent.${params.series}.gff"
    bedtools intersect -s -f 1.0 -r -loj -a ${params.series}.novel.gff -b  ${params.organism}.${params.release}.gff > novel_overlap_known_stringent.${params.series}.gff
    date
    fi

    ########################################
    
    if [ ! -s ${params.sajr_output}novel_overlap_known_stringent_novel.${params.series}.tsv ] ; then 
    # step 10
    # annotate_novel_segments.pl novel_overlap_known_stringent.gff > novel_overlap_known_stringent_novel.tsv
    
    cd ${params.sajr_code}
    echo "step 10: annotate novel segments"
    echo "perl annotate_novel_segments.pl ${params.sajr_output}novel_overlap_known_stringent.${params.series}.gff > ${params.sajr_output}novel_overlap_known_stringent_novel.${params.series}.tsv"
    perl annotate_novel_segments.pl ${params.sajr_output}novel_overlap_known_stringent.${params.series}.gff > ${params.sajr_output}novel_overlap_known_stringent_novel.${params.series}.tsv
    date
    fi

    ################################################################
    # The final part of the protocol is estimating inclusion levels
    # in each sample, and testing for differences between groups of
    # samples.
    ################################################################
    
    cd ${params.sajr_code}
    # step 11
    # Run SAJR in count mode for each sample using the novel.gff reference.
    
    echo "step 11: count reads with SAJR, this might take some time"
    echo "java -Xmx16G -jar sajr.jar count_reads ${params.tmp}/sajr.${params.series}.config"
    java -Xmx16G -jar sajr.jar count_reads ${params.tmp}/sajr.${params.series}.config

    echo "###############################################"
    echo "all done"
    date

        """

}


process sajr_diff_splicing {
    stageInMode 'symlink'
    stageOutMode 'move'

    script: 

    """
    #!/usr/bin/Rscript
    library(SAJR)
    library(biomaRt)
    library(openxlsx)
    library(dplyr)
    ensembl=useMart("ensembl")
    ensembl = useDataset("${params.spec}_gene_ensembl",mart=ensembl)
    gene_map=getBM(attributes=c("ensembl_gene_id", 'external_gene_name'), mart = ensembl, useCache = FALSE)
    run_sajr <- function(grp1, grp1_samples, grp2, grp2_samples, gene_map){
      
      significance=0.05
      ann.gff='${params.sajr_output}${params.series}.novel.gff'
      
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
      novelty<-read.delim("${params.sajr_output}novel_overlap_known_stringent_novel.${params.series}.tsv", header=F)
      novelty<-novelty[!duplicated(novelty[,1]),]
      colnames(novelty)=c("segment","Novelty")
      rownames(novelty)=novelty[,'segment']
      
      names_from_overlap<-read.delim("${params.sajr_output}novel2known_from_overlap.${params.series}.tsv", header=F)
      rownames(names_from_overlap)=names_from_overlap[,1]
      colnames(names_from_overlap)=c("sajrSegID","GeneName")
      head(names_from_overlap)
      
      names_from_comp<-read.delim("${params.sajr_output}sajr.${params.series}.novel2known.tsv", header=F)
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

    if (!file.exists('${params.sajr_output}count_files/')) {
    dir.create('${params.sajr_output}count_files/')
    }
    setwd('${params.sajr_output}count_files/')
    samples = read.xlsx("/workdir/nf-star-test/sample_sheet.xlsx")
    groups = unique(samples[,2])
    for(g1 in 1:length(groups)){
      for(g2 in 1:length(groups)){
        if (g1 < g2) {
          print(paste('testing', groups[g1], 'and', groups[g2]))
          run_sajr(groups[g1], samples[samples[,2] == groups[g1],'sample_name'],  groups[g2],  samples[samples[,2] == groups[g2],'sample_name'], gene_map)
        }
      }
    }
    sink('${params.sajr_output}count_files/sessionInfo.txt' )
    sessionInfo()
    sink()

    """

}

process upload_paths {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
  """
    cd ${params.scripts}
    rm -rf upload.txt
  
    cd ${params.fastqc_output}
  
    for file in *.html; do echo "${params.fastqc_output}\${file}" >> ${params.scripts}upload.txt_; done
  
    echo "multiqc ${params.multiqcOut}multiqc_report.html" >> ${params.scripts}upload.txt_ 
  
    cd ${params.bw_output}
    for file in *.bw ; do echo "${params.bw_output}\${file}" >>  ${params.scripts}upload.txt_ ; done
  
    cd ${params.sajr_output}count_files/
    for file in *.results.xlsx ; do echo "${params.sajr_output}\${file})" >> ${params.scripts}upload.txt_ ; done
  
    cd ${params.scripts}
    uniq upload.txt_ upload.txt 
    rm upload.txt_


  """
}



workflow images {

  main:
  get_images()
}

workflow config_template {
  fill_config_template()
}

workflow sajr_processing {
  sajr()
}

workflow sajr_diff {
  sajr_diff_splicing ()
}

workflow upload {
  main:
    upload_paths()
}
