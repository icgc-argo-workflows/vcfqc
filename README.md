## vcfqc

**icgc-argo-workflows/vcfqc** is a reproducible bioinformatics workflow that can be used to obtain QC metrics from variant calls in VCF/BCF format. It has been created to support quality control efforts within ICGC-ARGO project. The aggregated QC metrics are formed to align with the [GA4GH WGS_Quality_Control_Standards](https://www.ga4gh.org/product/wgs-quality-control-standards/).  

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

The workflow is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. 
The workflow has adopted [nf-core](https://nf-co.re/) framework and best practice guidelines to ensure reproducibility, portability and scalability. Where possible, many processes have been installed from [nf-core/modules](https://github.com/nf-core/modules). Moreover, ICGC ARGO specific modules have been installed form [icgc-argo-workflows/argo-modules](https://github.com/icgc-argo-workflows/argo-modules), which hosts ARGO reusable modules across all ICGC ARGO pipelines!

## Requirements

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Stage the required [reference files](#references)

## Quick start
1. Test the workflow running in `Local` mode on a minimal dataset with a single command:

   ```bash
   nextflow run icgc-argo-workflows/vcfqc \
     -profile test,standard \
     --outdir <OUTDIR>
   ```

2. Test the workflow running in `RDPC` mode with a single command if you have access to `RDPC-QA` env and have your valid api_token available:

   ```bash
   nextflow run icgc-argo-workflows/vcfqc \
     -profile test_rdpc_qa,standard,rdpc_qa \
     --api_token <YOUR_API_TOKEN> \
     --reference_base <REFERENCE_BASE> \
     --outdir <OUTDIR>
   ```

## Usage

### Workflow summary
Depending on where the input data are coming from and output data are sending to, the workflow can be running in two modes: `Local` and `RDPC` . The major tasks performed in the workflow are:
- (`RDPC` mode only) Download input sequencing metadata/data from data center using SONG/SCORE client tools
- Perform `Bcftools view` to count indels and snps counts and ratios
- Perform `Bcftools stats`  to collect stats for VCF
- (`RDPC` mode only) Generate SONG metadata for all collected QC metrics files and upload QC files to SONG/SCORE


### References
- Reference genome: 
  - GRCh38 reference genome fasta file. The file can be downloaded by:
    ```bash
    wget https://object.genomeinformatics.org/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa
    ``` 

  - GRCh38 reference genome fasta index file. The file can be downloaded by:
    ```bash
    wget https://object.genomeinformatics.org/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa.fai
    ```

- Autosome non-gap regions
  - `autosome_non_gap` bed file was downloaded from [NPM-sample-qc](https://raw.githubusercontent.com/c-BIG/NPM-sample-qc/master/resources/autosomes_non_gap_regions.bed) and staged under project folder [assets](https://github.com/icgc-argo-workflows/dnaalnqc/tree/main/assets)

> **NOTE**
> Please stage the reference files into the reference directory <REFERENCE_BASE> with the following folder structure
```bash
<REFERENCE_BASE>
├── GRCh38_hla_decoy_ebv.fa
├── GRCh38_hla_decoy_ebv.fa.fai
```


### Inputs
#### Local mode
First, prepare a sample sheet with your input data that looks as following example:

`sample_sheet.csv`:

```csv
sample,vcf
HG00100,assets/test/HG00100.hard-filtered.vcf.gz
HG00844,assets/test/HG00844.hard-filtered.vcf.gz
HG03722,assets/test/HG03722.hard-filtered.vcf.gz
```

Each row represents an aligned VCF or BCF from a sample.

Then, you need to download [Autosome non-gap regions](#references), and optionally [reference files](#references) staged in <REFERENCE_BASE>

Now, you can run the workflow using:

```bash
nextflow run icgc-argo-workflows/vcfqc \
   -profile <standard/singularity> \
   --local_mode true \
   --input sample_sheet.csv
   --outdir <OUTDIR>
```

#### RDPC mode
You can run the workflow in RDPC mode by using:
```bash
nextflow run icgc-argo-workflows/vcfqc \
    -profile <rdpc,rdpc_qa,rdpc_dev>,<standard/singularity> \
    --local_mode false \
    --study_id <STUDY_ID> \
    --analysis_ids <ANALYSIS_IDS> \
    --api_token <YOUR_API_TOKEN>
    --outdir <OUTDIR>
```
#### With additional arguements
You can run the workflow in RDPC mode by using:
```bash
nextflow run icgc-argo-workflows/vcfqc \
   -profile <standard/singularity> \
   --local_mode true \
   --input sample_sheet.csv \
   --fasta <REFERENCE_BASE>/GRCh38_hla_decoy_ebv.fa \
   --fasta_fai <REFERENCE_BASE>/GRCh38_hla_decoy_ebv.fa \
   --regions autosomes_non_gap_regions.bed \
   --outdir <OUTDIR>
```

> **NOTE**
> Please provide workflow parameters via the CLI or Nextflow `-params-file` option. 

### Outputs
Upon completion, you can find the aggregated QC metrics under directory:
```
/path/to/outdir/prep_metrics/<sample_id>.vcfqc.argo_metrics.json
/path/to/outdir/prep_metrics/<sample_id>.vcfqc.metrics.json
```

## Credits

icgc-argo-workflows/vcfqc was mostly written by Edmund Su (@edsu7), with contributions from 
Andrej Benjak, Charlotte Ng, Desiree Schnidrig, Linda Xiang, Miguel Vazquez, Morgan Taschuk, Raquel Manzano Garcia, Romina Royo and ICGC-ARGO Quality Control Working Group.  

Authors (alphabetical)
- Andrej Benjak
- Charlotte Ng
- Desiree Schnidrig
- Edmund Su
- Linda Xiang
- Miguel Vazquez
- Morgan Taschuk
- Raquel Manzano Garcia
- Romina Royo

## Citations
<!-- If you use  icgc-argo-workflows/dnaalnqc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
