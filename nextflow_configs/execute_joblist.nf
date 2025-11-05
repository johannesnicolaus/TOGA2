#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.joblist = 'NONE'  // file containing jobs

if (params.joblist == "NONE"){
    println("Usage: nextflow execute_joblist.nf  --joblist [joblist file] -c [config file]")
    System.exit(2);
}

lines = Channel.fromPath(params.joblist).splitText()

process execute_jobs {

    errorStrategy 'retry'
    maxRetries {}

    input:
    val line

    // one line represents an independent command
    script:
    """
    ${line}
    """
}

workflow {
    execute_jobs(lines)
}