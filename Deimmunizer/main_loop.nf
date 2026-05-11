#!/usr/bin/env nextflow

binders_in = binders_in = Channel
    .fromPath("${params.run_name}/alive_binders/*/*.pdb")
    .map { pdb ->
        tuple(
            pdb,
            pdb.parent.resolve('metrics.csv'),
            pdb.parent.resolve('info.txt'))
    }
params.allotypes = file("${params.run_name}/input/data/HLA_alleles.txt").text.trim()
params.data = file("data")

process LoadBinders {
    fair true

    input:
    tuple path(pdb_path), path(metrics), path(info)

    output:
    path "${pdb_path.simpleName}", emit: binder_dir
    path "${pdb_path.simpleName}/${pdb_path.name}", emit: full_pdb
    path "${pdb_path.simpleName}/${pdb_path.simpleName}_binder.pdb", emit: binder_pdb
    path "${pdb_path.simpleName}/${pdb_path.simpleName}_binder.fasta", emit: binder_fasta
    path "${pdb_path.simpleName}/${pdb_path.simpleName}_target.fasta", emit: target_fasta
    path "${pdb_path.simpleName}/metrics.csv", emit: metrics
    path "${pdb_path.simpleName}/info.txt", emit: info

    script:
    """
    mkdir ${pdb_path.simpleName}
    cp ${pdb_path} ${pdb_path.simpleName}/${pdb_path}
    echo ${pdb_path} | python ../../../modules/01loading.py
    cp "${metrics}" "./${pdb_path.simpleName}/metrics.csv"
    cp "${info}" "./${pdb_path.simpleName}/info.txt"
    """
}



process Mutate_py {
    fair true

    input:
    path pdb_path
    path binder_fasta
    path metrics_csv
    path target_fasta
    path info_parent
    path data

    output:
    path "./${pdb_path.simpleName}I/${pdb_path.simpleName}I_chain_A", emit: new_binder_fasta
    path "./${pdb_path.simpleName}I/info.txt", emit: info

    script:
    """
    mkdir ${pdb_path.simpleName}I
    cp "${info_parent}" "./${pdb_path.simpleName}I/info.txt"
    echo "${params.run_name}" | python ../../../modules/04mutate.py -m ${metrics_csv} -f ${binder_fasta} -t ${target_fasta} >> ${pdb_path.simpleName}I/info.txt
    """
}

process Binder_scoring {
    publishDir "${params.run_name}/new_gen", mode: 'copy'
    fair true    
    maxForks 1

    input:
    path new_binder
    path target_fasta
    path pdb_path //only for naming
    path info
    path metrics
    path data

    output:
    path "${pdb_path.simpleName}I/"
    path "${pdb_path.simpleName}I/${pdb_path.simpleName}I.pdb", emit: new_pdb

    script:
    """
    flock ../../../"${params.run_name}"/tmp/nextflow_gpu.lock -c '
        mkdir ${pdb_path.simpleName}I
        mkdir outDir
        echo "${params.run_name}" | python ../../../modules/02scoring.py -t ${target_fasta} -a ${new_binder}        
        
        (
            metrics=scoring_input.csv
            metrics_csv="\$metrics" ../../../modules/binder_scoring_seq_only.sh
        )
    '
    mv outDir/ipsae_and_ipae.csv ${pdb_path.simpleName}I/ipsae_and_ipae.csv 
    mv outDir/AF3/pdbs/*.pdb ${pdb_path.simpleName}I/${pdb_path.simpleName}I.pdb
    mv ${info} ${pdb_path.simpleName}I/${info}
    mv "${metrics}" ${pdb_path.simpleName}I/metrics.csv
   """
}

process BepiPred {
    fair true
    maxForks 1

    input:
    path input_fasta


    output:
    path "BepiPredOut/raw_output.csv"

    script:
    """
    flock ../../../"${params.run_name}"/tmp/nextflow_gpu.lock -c '
        mkdir BepiPredOut
        (
            source /dtu/projects/RFdiffusion/closed-loop/BepiPred3_src/Bepipred3_env/bin/activate
            python /dtu/projects/RFdiffusion/closed-loop/BepiPred3_src/bepipred3_CLI.py -i ${input_fasta} -o BepiPredOut -esm_dir ./ -pred vt_pred
        )
    '
    """
}

process DiscoTope {
    fair true

    input:
    path new_pdb

    output:
    path "./disco_out/output/${new_pdb.simpleName}_A_discotope3.csv"

    script:
    """
    mkdir disco_out
    (
        ../../../modules/DiscoTope.sh ${new_pdb}
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
    publishDir "${params.run_name}/new_gen", mode: 'copy'

    input:
    path new_pdb
    path binder_fasta
    path bepipred
    path netmhci
    path mhci_other
    path disco
    path data

    output:
    path "${new_pdb.simpleName}/metrics.csv", emit: full_metrics
    path "${new_pdb.simpleName}/metrics_epitopes.csv", emit: epitope_metrics


    script:
    """
    mkdir -p "${new_pdb.simpleName}"
    echo "${params.run_name}" | python ../../../modules/03in_silico_metrics.py -d ${disco} -p ${new_pdb} -b ${bepipred} -n ${netmhci} -f ${binder_fasta} "${new_pdb.simpleName}/metrics.csv"
    """
}
// Binder_scoring(LoadBinders.out.binder_fasta, LoadBinders.out.target_fasta, LoadBinders.out.full_pdb)




workflow {
    LoadBinders(binders_in)
    Mutate_py(LoadBinders.out.full_pdb, LoadBinders.out.binder_fasta, LoadBinders.out.metrics, LoadBinders.out.target_fasta, LoadBinders.out.info, params.data)
    Binder_scoring(Mutate_py.out.new_binder_fasta, LoadBinders.out.target_fasta, LoadBinders.out.full_pdb, Mutate_py.out.info, LoadBinders.out.metrics, params.data)
    NetMHCI(Mutate_py.out.new_binder_fasta,params.allotypes)
    BepiPred(Mutate_py.out.new_binder_fasta)
    DiscoTope(Binder_scoring.out.new_pdb)
    In_silico_metrics_eval(Binder_scoring.out.new_pdb ,Mutate_py.out.new_binder_fasta, BepiPred.out, NetMHCI.out.MHCI, NetMHCI.out.other, DiscoTope.out, params.data) 
    }