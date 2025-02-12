/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
// WorkflowPrealnqc.initialise(params, log)

// Check input path parameters to see if they exist
// def checkPathParamList = [ params.input, params.multiqc_config]

// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { INPUT_CHECK       } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
*/

// Set Channels for bcftools view
intervals_view_ch = params.regions ? Channel.fromPath(params.regions,checkIfExists: true).collect() : Channel.fromPath(params.autosome_non_gap,checkIfExists: true).collect()
targets_view_ch = params.targets ? Channel.fromPath(params.targets,checkIfExists: true).collect() : Channel.value([])
samples_view_ch = params.samples ? Channel.fromPath(params.samples,checkIfExists: true).collect() : Channel.value([])

// Set Channels for bcftools stats
intervals_stats_ch = params.regions ? Channel.fromPath(params.regions,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]} : Channel.fromPath(params.autosome_non_gap,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]}
targets_stats_ch = params.targets ? Channel.fromPath(params.targets,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]} : Channel.value([{},[]])
samples_stats_ch = params.samples ? Channel.fromPath(params.samples,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]} : Channel.value([{},[]])
exon_stats_ch = params.exons ? Channel.fromPath(params.exons,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]} : Channel.value([{},[]])
fasta_stats_ch = params.fasta ? Channel.fromPath(params.fasta,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]} : Channel.value([{},[]])
fasta_fai_stats_ch = params.fasta_fai ? Channel.fromPath(params.fasta_fai,checkIfExists: true).map{filePath -> [filePath.baseName,file(filePath)]} : Channel.value([{},[]])


/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT ARGO MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main' 
include { STAGE_INPUT       } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { CLEANUP           } from '../modules/icgc-argo-workflows/cleanup/main'
include { PREP_METRICS      } from '../modules/icgc-argo-workflows/prep/metrics/main'
include { PAYLOAD_QCMETRICS } from '../modules/icgc-argo-workflows/payload/qcmetrics/main'

include { BCFTOOLS_VIEW as BCF_SNP } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCF_INS } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCF_DEL } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCF_HET_SNP } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCF_HOMO_SNP } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCF_HET_INDEL } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCF_HOMO_INDEL } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_STATS as STATS_TSTV } from '../modules/nf-core/bcftools/stats/main'

//include { BCFTOOLS_VIEW as BCF_TI } from '../modules/nf-core/bcftools/view/main'
//include { BCFTOOLS_VIEW as BCF_TA } from '../modules/nf-core/bcftools/view/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCFQC {

    ch_versions = Channel.empty()

    // Stage input files
    ////Integration with song works but put on hold for now until support for local is better
    STAGE_INPUT(params.study_id, params.analysis_ids, params.input)
    ch_versions = ch_versions.mix(STAGE_INPUT.out.versions)


    //https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#count-snvs
    ////Count PASS SNPs
    BCF_SNP(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )
    ch_versions = ch_versions.mix(BCF_SNP.out.versions)
    BCF_SNP_COUNT =  BCF_SNP.out.vcf.map{
        meta,vcf -> 
                vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.snp.vcf")
        return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.snp.vcf")]
        }

    ////https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#count-insertions
    ////Count PASS insertions
    BCF_INS(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )

    ch_versions = ch_versions.mix(BCF_INS.out.versions)

    BCF_INS_COUNT =  BCF_INS.out.vcf.map{
    meta,vcf -> 
            vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.ins.vcf")
    return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.ins.vcf")]
    }
    
    ////https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#count-deletions
    ////Count Pass deletions
    BCF_DEL(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )

    ch_versions = ch_versions.mix(BCF_DEL.out.versions)

    BCF_DEL_COUNT =  BCF_DEL.out.vcf.map{
    meta,vcf -> 
            vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.del.vcf")
    return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.del.vcf")]
    }


    ////https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#ratio-insertionsdeletions
    ////Ratio of PASS insertsions/PASS deletions
    //Calculated in PREP_METRICS using BCF_INS_COUNT and BCF_DEL_COUNT

    ///https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#ratio-heterozygoushomozygous-snvs
    ////Count PASS heterozygous SNPs
    BCF_HET_SNP(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )
    ch_versions = ch_versions.mix(BCF_HET_SNP.out.versions)
    BCF_HET_SNP_COUNT =  BCF_HET_SNP.out.vcf.map{
        meta,vcf -> 
                vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.het_snp.vcf")
        return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.het_snp.vcf")]
        }

    ////Count PASS homozygous SNPs
    BCF_HOMO_SNP(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )
    ch_versions = ch_versions.mix(BCF_HOMO_SNP.out.versions)
    BCF_HOMO_SNP_COUNT =  BCF_HOMO_SNP.out.vcf.map{
        meta,vcf -> 
                vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.homo_snp.vcf")
        return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.homo_snp.vcf")]
        }

    ////Ratio of PASS heterozygous SNPs/PASS homozygous SNPs
    //Calculated in PREP_METRICS using BCF_HET_SNP_COUNT and BCF_HOMO_SNP_COUNT

    ////https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#ratio-heterozygoushomozygous-indels
    ////Count PASS heterozygous indels
    BCF_HET_INDEL(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )
    ch_versions = ch_versions.mix(BCF_HET_INDEL.out.versions)
    BCF_HET_INDEL_COUNT =  BCF_HET_INDEL.out.vcf.map{
        meta,vcf -> 
                vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.het_indel.vcf")
        return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.het_indel.vcf")]
        }

    ////Count PASS homozygous indels
    BCF_HOMO_INDEL(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_view_ch , // regions
        targets_view_ch , //target
        samples_view_ch  //samples
    )
    ch_versions = ch_versions.mix(BCF_HOMO_INDEL.out.versions)
    BCF_HOMO_INDEL_COUNT =  BCF_HOMO_INDEL.out.vcf.map{
        meta,vcf -> 
                vcf.renameTo("${vcf.parent}/${meta.patient}.${meta.sample}.homo_indel.vcf")
        return [meta,file("${vcf.parent}/${meta.patient}.${meta.sample}.homo_indel.vcf")]
        }

    ////Ratio of PASS heterozygous indels/PASS homozygous indels
    //Calculated in PREP_METRICS using BCF_HET_INDEL_COUNT and BCF_HOMO_INDEL_COUNT

    ////https://github.com/ga4gh/quality-control-wgs/blob/main/metrics_definitions/metrics_definitions.md#ratio-transitionstransversions-titv
    ////Ratio of PASS transitions/PASS transversions
    STATS_TSTV(
        STAGE_INPUT.out.meta_files, // meta,vcf,tbi
        intervals_stats_ch , // meta,regions
        targets_stats_ch  , // meta,target
        samples_stats_ch , // meta,samples
        exon_stats_ch , // meta,exons
        fasta_stats_ch   // meta,fasta
    )
    ch_versions = ch_versions.mix(STATS_TSTV.out.versions)


    ////SANITY CHECK
    //BCF_SNP_COUNT.subscribe{println "SNP ${it}"}
    //BCF_INS_COUNT.subscribe{println "INS ${it}"}
    //BCF_DEL_COUNT.subscribe{println "DEL ${it}"}
    //RATIO_INS_DEL.subscribe{println "RATIO_INS_DEL ${it}"}
    //BCF_HET_SNP_COUNT.subscribe{println "HET_SNP ${it}"}
    //BCF_HOMO_SNP_COUNT.subscribe{println "HOMO_SNP ${it}"}
    //RATIO_HET_HOMO_SNP.subscribe{println "RATIO_HET_HOMO_SNP ${it}"}
    //BCF_HET_INDEL_COUNT.subscribe{println "HET_INDEL ${it}"}
    //BCF_HOMO_INDEL_COUNT.subscribe{println "HOMO_INDEL ${it}"}
    // RATIO_HET_HOMO_INDEL.subscribe{println "RATIO_HET_HOMO_INDEL ${it}"}
    // RATIO_TI_TA.subscribe{println "TI_TA ${it}"}

    //Collect VCF counts for Prep metrics
    ch_prep_metrics=STATS_TSTV.out.stats
    .join(BCF_SNP_COUNT)
    .join(BCF_INS_COUNT)
    .join(BCF_DEL_COUNT)
    .join(BCF_HET_SNP_COUNT)
    .join(BCF_HOMO_SNP_COUNT)
    .join(BCF_HET_INDEL_COUNT)
    .join(BCF_HOMO_INDEL_COUNT)
    .map{
        meta,fileA,fileB,fileC,fileD,fileE,fileF,fileG,fileH ->
        [meta,[fileA,fileB,fileC,fileD,fileE,fileF,fileG,fileH]]
    }

    PREP_METRICS(
        ch_prep_metrics.map{ meta,files -> [meta,[]]},
        ch_prep_metrics.map{ meta,files -> files}
    )

    // Combine channels to determine upload status and payload creation
    // make metadata and files match  
    STAGE_INPUT.out.meta_analysis.map { meta, metadata -> [[id: meta.sample, study_id: meta.study_id], metadata]}
        .unique().set{ ch_meta_metadata }
  
    ch_meta_metadata
    .join(STATS_TSTV.out.stats.map{ meta, stats -> [[id:meta.sample, study_id: meta.study_id],stats]})
    .join(PREP_METRICS.out.metrics_json.map{ meta, metrics -> [[id:meta.sample, study_id: meta.study_id],metrics]})
    .set { ch_metadata_files }


    STAGE_INPUT.out.upRdpc.combine(ch_metadata_files)
    .map{upRdpc, meta, metadata, files, metrics -> 
    [[id: meta.id, study_id: meta.study_id, upRdpc: upRdpc],
      files,metadata, metrics]}
    .branch{
      upload: it[0].upRdpc
    }.set{ch_metadata_files_status}

    //Collect Software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))


    //Make payload
    PAYLOAD_QCMETRICS(
        ch_metadata_files_status.upload, CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect()
        )
    //Upload payload
    SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)
    
    //cleanup the files if specified
    if (params.cleanup) {
      // Gather files to remove   
      ch_files = Channel.empty()
      ch_files = ch_files.mix(STAGE_INPUT.out.meta_files) 
      ch_files.map{ meta, file1, file2 -> [file1, file2]}
      .set { ch_files_to_remove1 }

      PAYLOAD_QCMETRICS.out.payload_files
      .map {meta, payload, files -> files}
      .unique()
      .set { ch_files_to_remove2 }

      ch_files_to_remove = Channel.empty()
      ch_files_to_remove = ch_files_to_remove.mix(BCF_SNP_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_INS_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_HET_SNP_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_HOMO_SNP_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_HET_SNP_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_HET_INDEL_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_HOMO_INDEL_COUNT)
      ch_files_to_remove = ch_files_to_remove.mix(BCF_HOMO_INDEL_COUNT)

      CLEANUP(ch_files_to_remove.unique().collect(), SONG_SCORE_UPLOAD.out.analysis_id)    
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/