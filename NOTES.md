## Development Notes

This document contains notes for development; goals, to do's, ect..

### Goals

The aim of this pipline is to centralize the development of pipelines used for all things viral deepsequencing. The key principle is that you shoul be able to download data from the SRA Run Selector and use that csv, with other added information, to run the whole pipeline. The user should be able to toggle the config file to determine which aligners, quality control, and analysis tools they want to use. 

### To Do

- I now have all of the varaint calling rules worked into the new framework that can use multiple viruses or hosts. I still need to make sure all of the qc is updated, but otherwise I can now start looking into new modules like assembly. 

- I was looking through some variant calling pipelines for RNA-seq and I was looking through Maddy's pipeline and I decided to try some different approachs to incorporating `GATK` best practices. I can't get Maddy's BSQR (Base Score Recallibration) rule to work for my pipeline. There is no 'know-sites' to use for viruses. I'm not sure it makes sense to do any kind of recallibration for this analysis? I can't find anywhere that uses realignment and indel qualities for small genomes with high error rates. 

- I was reading about some of the caveats of variant calling with viruses. I should look into incorporating 'DeepSNV'. I think it's an R-script and it was benchmarked on viruses against lofreq by J.T. Macrone. 

- `Lofreq` can produce an empty VCF file if it finds no variants. It doesn't get treated as an empty file by the pipeline because it had a header. It can cause the `vcf_to_table` rule to fail. Possible putting in checkpoints here would solve this issue, or splitting the shell command with try/except in python would fix the issue.

- Add in a rule to index genomes with samtools before variant calling. Lofreq can break if STAR alignment happens before BWA (because the genome used to call variants won't have been indexed yet).  

- With BWA at least, there are some reads pairs that map to different chromosomes. This proves to be an issue for Lofreq. It only seems to happen with BWA. It's very curious. 

### Tools to Use

Here are some of the functionalties that I'm trying to build into the pipeline and the tools that can be used to make this work.

- Taxonomy of input samples 

- 
