# nf-sajr


Run the workflow:

```
nextflow run nf-sajr -params-file nf-sajr/params.json -entry images -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry config_template -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry sajr_processing -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry sajr_diff -profile studio && \
nextflow run nf-sajr -params-file nf-sajr/params.json -entry upload -profile studio

___


## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag> 
```
___


