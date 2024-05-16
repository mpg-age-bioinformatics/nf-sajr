# nf-sajr

```
nextflow run nf-sajr -params-file nf-sajr/params.json -entry images -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry config_template -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry sajr_processing -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry sajr_diff -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry upload -profile studio

```
