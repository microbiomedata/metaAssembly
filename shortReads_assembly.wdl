version 1.0

struct Resources {
    Float kmers
    Int predicted
    Int request
    Int cpu
    Int runtime_minutes
}

workflow jgi_metaASM {
    input {
        # String? outdir
        String? input_file
        String proj
        String prefix=sub(proj, ":", "_")
        String rename_contig_prefix="scaffold"
        # Float uniquekmer=1000
        String bbtools_container="microbiomedata/bbtools:39.03"
        String spades_container="staphb/spades:4.0.0"
        String workflowmeta_container="microbiomedata/workflowmeta:1.1.1"
        Boolean paired = true
        Int     stage_mem = 4       # GB
        Int     stage_cpu = 2
        Int     stage_run_mins = 20
        Int     bbcms_cpu = 40
        Int     bbcms_mem = 120     # GB
        Int     bbcms_run_mins = 45
        String  pred_mem = "100MB"
        Int     pred_cpu = 1
        Int     pred_run_mins = 5
        Int     agp_mem = 10        # GB
        Int     agp_cpu = 1
        Int     agp_run_mins = 5
        Int     mapping_mem = 120   # GB
        Int     mapping_cpu = 32
        Int     mapping_run_mins = 60
        String  make_info_mem = "100MB"
        Int     make_info_cpu = 1
        Int     make_info_run_mins = 5
        String  finish_mem = "100MB"
        Int     finish_cpu = 2
        Int     finish_run_mins = 30
    }

    call stage {
        input:
            input_file=input_file,
            container=bbtools_container,
            memory=stage_mem,
            cpu = stage_cpu,
            run_mins = stage_run_mins
    }

    call bbcms {
        input: 
            input_files=stage.assembly_input, 
            paired=paired,
            container=bbtools_container, 
            memory= bbcms_mem,  
            cpu = bbcms_cpu,
            run_mins = bbcms_run_mins
    }

    call predict_memory {
        input:
            kmer_count = bbcms.floatkmer,
            container = bbtools_container,
            memory=pred_mem,
            cpu = pred_cpu,
            run_mins = pred_run_mins
    }

    call assy {
        input: 
            infile1=bbcms.out1, 
            infile2=bbcms.out2, 
            paired=paired,
            container=spades_container, 
            threads=predict_memory.resource.cpu, 
            memory=predict_memory.resource.request,
            cpu=predict_memory.resource.cpu, 
            run_mins = predict_memory.resource.runtime_minutes
    }

    call create_agp {
        input: 
            scaffolds_in=assy.out, 
            rename_contig_prefix = rename_contig_prefix, 
            container=bbtools_container, 
            memory=agp_mem,
            cpu = agp_cpu,
            run_mins = agp_run_mins
    }

    call read_mapping_pairs {
        input: 
            reads=stage.assembly_input, 
            ref=create_agp.outcontigs, 
            paired=paired,
            container=bbtools_container, 
            memory=mapping_mem, 
            cpu = mapping_cpu,
            run_mins = mapping_run_mins
    }

    call make_info_file {
        input: 
            bbcms_info=bbcms.outcounts, 
            assy_info=assy.outlog, 
            prefix=prefix,
            container=bbtools_container, 
            memory=make_info_mem,
            cpu = make_info_cpu,
            run_mins = make_info_run_mins
    }

    call finish_asm {
        input:
            proj=proj,
            prefix=prefix,
            # start=stage.start, 
            fasta=create_agp.outcontigs,
            scaffold=create_agp.outscaffolds,
            agp=create_agp.outagp,
            bam=read_mapping_pairs.outbamfile,
            samgz=read_mapping_pairs.outsamfile,
            covstats=read_mapping_pairs.outcovfile,
            asmstats=create_agp.outstats,
            bbcms_fastq = bbcms.out,
            container=workflowmeta_container,
            memory=finish_mem,
            cpu = finish_cpu,
            run_mins = finish_run_mins
    }
 
    output {
        File contig=finish_asm.outcontigs
        File scaffold=finish_asm.outscaffolds
        File agp=finish_asm.outagp
        File bam=finish_asm.outbam
        File samgz=finish_asm.outsamgz
        File covstats=finish_asm.outcovstats
        File asmstats=finish_asm.outasmstats
        File asminfo=make_info_file.asminfo
        File bbcms_fastq = finish_asm.outbbcms
    }
 
    meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.0"
    }

}

task stage {
   input { 
        String? input_file
        String target = "staged.fastq.gz"
        String output1 = "input.left.fastq.gz"
        String output2 = "input.right.fastq.gz"
        String container
        Int    memory
        Int    xmxmem = floor(memory * 0.75)
        Int    cpu
        Int    run_mins
   }

   command <<<
       set -euo pipefail
       if [ "$( echo ~{input_file} | grep -E -c 'https\?:' )" -gt 0 ] ; then
           wget ~{input_file} -O ~{target}
       else
           ln -s ~{input_file} ~{target} || cp ~{input_file} ~{target}
       fi

        reformat.sh \
        ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx10G" }\
        in=~{target} \
        out1=~{output1} \
        out2=~{output2}    
       # Capture the start time
       date --iso-8601=seconds > start.txt

   >>>

   output {
      Array[File] assembly_input = [output1, output2]
      String start = read_string("start.txt")
   }
   runtime {
        docker: container
        memory: memory
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
   }
}

task make_info_file {
    input{
        File assy_info
        File bbcms_info
        String prefix
        String container
        String memory
        Int    cpu
        Int    run_mins
    }

    command<<<
        set -euo pipefail
        bbtools_version=$( grep BBToolsVer ~{bbcms_info} | awk '{print $2}' | sed -e 's/"//g' -e 's/,//' ) 
        spades_version=$( grep  'SPAdes version' ~{assy_info} | awk '{print $3}' )
        echo -e "The workflow takes paired-end reads runs error correction by bbcms.sh (BBTools(1) version $bbtools_version)." > ~{prefix}_metaAsm.info
        echo -e "The clean reads are assembled by metaSpades(2) version $spades_version with parameters, --only-assembler -k 33,55,77,99,127  --meta" >> ~{prefix}_metaAsm.info
        echo -e "After assembly, Contigs and Scaffolds are consumed by the *create_agp* task to rename the FASTA header and generate an AGP format (https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) file which describes the assembly"  >> ~{prefix}_metaAsm.info
        echo -e "In the end, the reads are mapped back to contigs by bbmap (BBTools(1) version $bbtools_version) for coverage information." >> ~{prefix}_metaAsm.info

        echo -e "\n(1) B. Bushnell: BBTools software package, http://bbtools.jgi.doe.gov/" >> ~{prefix}_metaAsm.info
        echo -e "(2) Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. Genome Res. 2017 May;27(5):824-834."  >> ~{prefix}_metaAsm.info
    >>>

    output {
        File asminfo = "~{prefix}_metaAsm.info"
    }
    runtime {
        docker: container
        memory: memory
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
    }
}

task finish_asm {
    input {
        File fasta
        File scaffold
        File? agp
        File bam
        File? samgz
        File? covstats
        File asmstats
        File bbcms_fastq
        String proj
        String prefix 
        String orig_prefix="scaffold"
        String sed="s/~{orig_prefix}_/~{proj}_/g"
        String container
        String memory
        Int    cpu
        Int    run_mins
    }

    command<<<
        set -euo pipefail
        # end=$( date --iso-8601=seconds )

        ## RE-ID
        cat ~{fasta} | sed ~{sed} > ~{prefix}_contigs.fna
        cat ~{scaffold} | sed ~{sed} > ~{prefix}_scaffolds.fna
        cat ~{covstats} | sed ~{sed} > ~{prefix}_covstats.txt
        cat ~{agp} | sed ~{sed} > ~{prefix}_assembly.agp
        ln ~{bbcms_fastq} ~{prefix}_bbcms.fastq.gz || ln -s ~{bbcms_fastq} ~{prefix}_bbcms.fastq.gz

        ## Bam file     
        samtools view -h ~{bam} | sed ~{sed} | \
          samtools view -hb -o ~{prefix}_pairedMapped_sorted.bam

        ## Sam.gz file
        samtools view -h ~{samgz} | sed ~{sed} | \
          gzip -c - > ~{prefix}_pairedMapped.sam.gz

        ## Remove an extra field from the stats and normalize keys
        cat ~{asmstats} | jq 'del(.filename)' > stats.json.tmp

        python3 <<EOF
        import json

        with open("stats.json.tmp") as f:
            stats = json.load(f)

        key_map = {
            "ctg_N50": "ctg_n50",
            "ctg_L50": "ctg_l50",
            "ctg_N90": "ctg_n90",
            "ctg_L90": "ctg_l90",
            "scaf_N50": "scaf_n50",
            "scaf_L50": "scaf_l50",
            "scaf_N90": "scaf_n90",
            "scaf_L90": "scaf_l90",
            "scaf_n_gt50K": "scaf_n_gt50k",
            "scaf_l_gt50K": "scaf_l_gt50k",
            "scaf_pct_gt50K": "scaf_pct_gt50k"
        }

        for old, new in key_map.items():
            if old in stats:
                stats[new] = stats.pop(old)

        with open("stats.json", "w") as f:
            json.dump(stats, f, indent=2)
        EOF
    >>>

    output {
        File outcontigs = "~{prefix}_contigs.fna"
        File outscaffolds = "~{prefix}_scaffolds.fna"
        File outagp = "~{prefix}_assembly.agp"
        File outbam = "~{prefix}_pairedMapped_sorted.bam"
        File outsamgz = "~{prefix}_pairedMapped.sam.gz"
        File outcovstats = "~{prefix}_covstats.txt"
        File outasmstats = "stats.json"
        File outbbcms = "~{prefix}_bbcms.fastq.gz"
    }

    runtime {
        docker: container
        memory: memory
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
    }
}

task read_mapping_pairs{
    input {
    Array[File] reads
    File ref
    String container
    Int     memory
    Int     xmxmem = floor(memory * 0.75)
    Int     cpu
    Int     threads = cpu * 2
    Int     run_mins
    Boolean paired = true
    String bbmap_interleaved_flag = if paired then 'interleaved=true' else 'interleaved=false'

    String filename_unsorted="pairedMapped.bam"
    String filename_outsam="pairedMapped.sam.gz"
    String filename_sorted="pairedMapped_sorted.bam"
    String filename_sorted_idx="pairedMapped_sorted.bam.bai"
    String filename_bamscript="to_bam.sh"
    String filename_cov="covstats.txt"
    String system_cpu="$( grep \"model name\" /proc/cpuinfo | wc -l )"
    String jvm_threads=select_first([threads,system_cpu])
    }

    command<<<
        set -euo pipefail
        if [[ ~{reads[0]}  == *.gz ]] ; then
             cat ~{sep=" " reads} > infile.fastq.gz
             export mapping_input="infile.fastq.gz"
        fi
        if [[ ~{reads[0]}  == *.fastq ]] ; then
             cat ~{sep=" " reads} > infile.fastq
             export mapping_input="infile.fastq"
        fi

        bbmap.sh \
        ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx105G" } \
        threads=~{jvm_threads} \
        nodisk=true \
        ~{bbmap_interleaved_flag} \
        ambiguous=random \
        in="$mapping_input" \
        ref=~{ref} \
        out=~{filename_unsorted} \
        covstats=~{filename_cov} \
        bamscript=~{filename_bamscript}

        samtools sort \
        -m100M \
        -@ \
        ~{jvm_threads} \
        ~{filename_unsorted} \
        -o ~{filename_sorted}

        samtools index ~{filename_sorted}

        reformat.sh \
        ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx105G" } \
        in=~{filename_unsorted} \
        out=~{filename_outsam} \
        overwrite=true

        ln -s ~{filename_cov} mapping_stats.txt
        rm "$mapping_input"

    >>>

    output {
        File outbamfile = filename_sorted
        File outbamfileidx = filename_sorted_idx
        File outcovfile = filename_cov
        File outsamfile = filename_outsam
    }

    runtime {
        docker: container
        memory: memory
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
    }
}

task create_agp {
    input {
    File scaffolds_in
    String rename_contig_prefix
    String prefix="assembly"
    String filename_contigs="~{prefix}_contigs.fna"
    String filename_scaffolds="~{prefix}_scaffolds.fna"
    String filename_agp="~{prefix}.agp"
    String filename_legend="~{prefix}_scaffolds.legend"
    String container
    Int    memory
    Int    xmxmem = floor(memory * 0.75)
    Int    cpu
    Int    run_mins
    }

    command<<<
        set -euo pipefail
        fungalrelease.sh \
        ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx105G" } \
        in=~{scaffolds_in} \
        out=~{filename_scaffolds} \
        outc=~{filename_contigs} \
        agp=~{filename_agp} \
        legend=~{filename_legend} \
        mincontig=200 \
        minscaf=200 \
        sortscaffolds=t \
        sortcontigs=t \
        overwrite=t

        if [ "~{rename_contig_prefix}" != "scaffold" ]; then
            sed -i 's/scaffold/~{rename_contig_prefix}_scf/g' \
            ~{filename_contigs} ~{filename_scaffolds} ~{filename_agp} ~{filename_legend}
        fi
        bbstats.sh format=8 in=~{filename_scaffolds} out=stats.json
        sed -i 's/l_gt50k/l_gt50K/g' stats.json

    >>>

    output {
    File outcontigs = filename_contigs
    File outscaffolds = filename_scaffolds
    File outagp = filename_agp
    File outstats = "stats.json"
    File outlegend = filename_legend
    }

    runtime {
        docker: container
        memory: "~{memory} GiB"
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
    }
}

task assy {
    input{
    File infile1
    File infile2
    String container
    Int? threads
    Int? memory
    Int? cpu
    Int run_mins
    String outprefix="spades3"
    String filename_outfile="~{outprefix}/scaffolds.fasta"
    String filename_spadeslog ="~{outprefix}/spades.log"
    String system_cpu="$( grep \"model name\" /proc/cpuinfo | wc -l )"
    String spades_cpu=select_first([threads,system_cpu])
    Boolean paired = true
    }

     command <<<
        set -euo pipefail
        if ~{paired}; then
            spades.py \
            -m 2000 \
            -o ~{outprefix} \
            --only-assembler \
            -k 33,55,77,99,127  \
            --meta \
            -t ~{spades_cpu} \
            -1 ~{infile1} \
            -2 ~{infile2}
        else
            spades.py \
            -m 2000 \
            -o ~{outprefix} \
            --only-assembler \
            -k 33,55,77,99,127 \
            -t ~{spades_cpu} \
            -s ~{infile1}
        fi
    >>>

    output {
        File out = filename_outfile
        File outlog = filename_spadeslog
    }

    runtime {
        docker: container
        memory: "~{memory} GiB"
        cpu: cpu
        runtime_minutes: run_mins
    }
}

task bbcms {
    input{
    Array[File] input_files
    String container
    Int    memory
    Int    xmxmem = floor(memory * 0.75)
    Int    cpu
    Int    run_mins
    Boolean paired = true
    String filename_outfile="input.corr.fastq.gz"
    String filename_outfile1="input.corr.left.fastq.gz"
    String filename_outfile2="input.corr.right.fastq.gz"
    String filename_readlen="readlen.txt"
    String filename_outlog="stdout.log"
    String filename_errlog="stderr.log"
    String filename_kmerfile="unique31mer.txt"
    String filename_counts="counts.metadata.json"
    }

    command<<<
        set -euo pipefail
        if file --mime -b ~{input_files[0]} | grep gzip; then
            cat ~{sep=" " input_files} > infile.fastq.gz
            export bbcms_input="infile.fastq.gz"
        fi

        if file --mime -b ~{input_files[0]} | grep plain; then
            cat ~{sep=" " input_files} > infile.fastq
            export bbcms_input="infile.fastq"
        fi

        bbcms.sh \
        ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx105G" } \
        metadatafile=~{filename_counts} \
        mincount=2 \
        highcountfraction=0.6 \
        in="$bbcms_input" \
        out=~{filename_outfile} \
        > >(tee -a ~{filename_outlog}) \
        2> >(tee -a ~{filename_errlog} >&2) \
        && grep Unique ~{filename_errlog} \
        | rev |  cut -f 1 | rev  \
        > ~{filename_kmerfile}
        
        if ~{paired}; then
            reformat.sh \
            ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx105G" } \
            in=~{filename_outfile} \
            out1=~{filename_outfile1} \
            out2=~{filename_outfile2}
        fi

        readlength.sh \
        ~{if (defined(memory)) then "-Xmx" + xmxmem + "G" else "-Xmx105G" } \
        in=~{filename_outfile} \
        out=~{filename_readlen}

        rm "$bbcms_input"
        
    >>>

    output {
        File out = filename_outfile
        File out1 = if paired then filename_outfile1 else filename_outfile
        File out2 = if paired then filename_outfile2 else filename_outfile
        File outreadlen = filename_readlen
        File stdout = filename_outlog
        File stderr = filename_errlog
        File outcounts = filename_counts
        File outkmer = filename_kmerfile
        Float floatkmer = read_float(filename_kmerfile)  
    }

    runtime {
        docker: container
        memory: "~{memory} GiB"
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
    }
}

task predict_memory {
    input {
        Float  kmer_count
        String json_out = "outfile.json"
        String container
        String memory
        Int    cpu
        Int    run_mins
    }

    command <<<
        python <<CODE
        import json
        import math
        kmers = ~{kmer_count}
        predicted_mem =  (kmers * 1.416e-08 + 8.676e-01) * 1.1
        predicted_time = (kmers * 2.153e-09 - 6.437e-01) * 1.5
        rounded_time = int(math.ceil(predicted_time / 10.0) * 10) * 60
        rounded_time = 10 if rounded_time <= 1 else rounded_time
        (mem, cpu)  = next(((m, c) for p, m, c in [(120, 110, 16), (250, 240, 32), (500, 490, 32)] if predicted_mem < p), (490, 32))
        with open("~{json_out}", "w") as f:
            f.write(json.dumps({"kmers": float(round(kmers)), "predicted": int(round(predicted_mem)), "request": int(mem), "cpu": int(cpu), "runtime_minutes": int(rounded_time)}))
        CODE
    >>>
    runtime {
        docker: container
        memory: memory
        cpu:  cpu
        runtime_minutes: run_mins
        maxRetries: 1
    }
    output {
        Resources resource = read_json(json_out)
    }
}