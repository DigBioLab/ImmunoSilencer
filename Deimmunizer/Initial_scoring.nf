#!/usr/bin/env nextflow
binders_in = channel.fromPath("${params.run_name}/input/*.pdb") //main.sh makes sure this is only one file
params.allotypes = file("${params.run_name}/input/data/HLA_alleles.txt").text.trim()
params.data = file("data")

process LoadBinders {

    input:
    path pdb_path

    output:
    path "${pdb_path.simpleName}", emit: binder_dir
    path "${pdb_path.simpleName}/${pdb_path.name}", emit: full_pdb
    path "${pdb_path.simpleName}/${pdb_path.simpleName}_binder.pdb", emit: binder_pdb
    path "${pdb_path.simpleName}/${pdb_path.simpleName}_binder.fasta", emit: binder_fasta
    path "${pdb_path.simpleName}/${pdb_path.simpleName}_target.fasta", emit: target_fasta

    script:
    """
    mkdir ${pdb_path.simpleName}
    cat ${pdb_path} > ${pdb_path.simpleName}/${pdb_path.name}
    echo ${pdb_path} | python ../../../modules/01loading.py
    """
}

process BepiPred {
    fair true
    maxForks 1

    input:
    path input_fasta
    path input_pdb

    output:
    path "BepiPredOut/raw_output.csv"

    script:
    """
    flock "../../../${params.run_name}"/tmp/nextflow_gpu.lock -c '
        mkdir "${params.run_name}/BepiPredOut"
        (
            source /dtu/projects/RFdiffusion/closed-loop/BepiPred3_src/Bepipred3_env/bin/activate
            python /dtu/projects/RFdiffusion/closed-loop/BepiPred3_src/bepipred3_CLI.py -i ${input_fasta} -o BepiPredOut -esm_dir ./ -pred vt_pred
        )
    '
    """
}

process DiscoTope{
    fair true

    input:
    path input_pdb

    output:
    path "disco_out/output/${input_pdb.simpleName}_A_discotope3.csv"

    script:
    """
    mkdir disco_out
    (
        ../../../modules/DiscoTope.sh ${input_pdb}
    )
    """
}

process NetMHCI {
    fair true
    input:
    path input_fasta
    val allotypes
    
    output:
    path "NetMHCI.out", emit: MHCI
    path "other.dat", emit: other // contains 3 lines: filename of binder.fasta, allotypes, "NetMHCI.out"  (name of outputfile of NetMHCI)


    script:
    """
    echo ${input_fasta} > other.dat
    echo ${allotypes} >> other.dat
    /dtu/projects/RFdiffusion/closed-loop/netMHCpan-4.2/netMHCpan -f ${input_fasta} -a ${allotypes} > NetMHCI.out
    echo NetMHCI.out >> other.dat
    """
}

process In_silico_metrics_eval {
    fair true


    input:
    path binder_dir
    path full_pdb
    path _binder_fasta
    path bepipred
    path netmhci
    path mhci_other
    path disco
    path data

    output:
    path "${binder_dir}/metrics.csv", emit: full_metrics
    path "${binder_dir}/metrics_epitopes.csv", emit: epitope_metrics


    script:
    """
    echo "${params.run_name}" | python ../../../modules/03in_silico_metrics.py -d ${disco} -p ${full_pdb} -b ${bepipred} -n ${netmhci} -f ${_binder_fasta} ${binder_dir}/metrics.csv > debug.txt
    """
}


process Binder_scoring {

    publishDir "${params.run_name}/input/data", mode: 'copy'
    
    input:
    path binder_fasta
    path target_fasta
    path pdb_path
    path metrics
    path epitope_metrics
    path data

    output:
    path "${pdb_path.simpleName}/"

    script:
    """
    mkdir ${pdb_path.simpleName}
    mkdir outDir
    echo "${params.run_name}" | python ../../../modules/02scoring.py -t ${target_fasta} -a ${binder_fasta} > debug.txt
      
    (
        metrics="scoring_input.csv"
        metrics_csv="\$metrics" ../../../modules/binder_scoring_seq_only.sh
    )

    mv outDir/ipsae_and_ipae.csv ${pdb_path.simpleName}/ipsae_and_ipae.csv 
    mv "${pdb_path}" "${pdb_path.simpleName}"/${pdb_path.name}
    touch ${pdb_path.simpleName}/info.txt
    mv "${epitope_metrics}" ${pdb_path.simpleName}/metrics_epitopes.csv
    mv "${metrics}" ${pdb_path.simpleName}/metrics.csv

    """
}
workflow {
    LoadBinders(binders_in)
    NetMHCI(LoadBinders.out.binder_fasta,params.allotypes)
    BepiPred(LoadBinders.out.binder_fasta, LoadBinders.out.binder_pdb)
    DiscoTope(LoadBinders.out.binder_pdb)
    In_silico_metrics_eval(LoadBinders.out.binder_dir,LoadBinders.out.full_pdb ,LoadBinders.out.binder_fasta, BepiPred.out, NetMHCI.out.MHCI, NetMHCI.out.other, DiscoTope.out,params.data ) 
    Binder_scoring(LoadBinders.out.binder_fasta, LoadBinders.out.target_fasta, LoadBinders.out.full_pdb, In_silico_metrics_eval.out.full_metrics, In_silico_metrics_eval.out.epitope_metrics, params.data)
}