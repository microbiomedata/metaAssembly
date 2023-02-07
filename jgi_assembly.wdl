workflow jgi_metaASM {
    String? memory
    String? threads
    String input_file
    String proj
    String rename_contig_prefix = "scaffold"
    Float uniquekmer = 1000
    String bbtools_container = "microbiomedata/bbtools:38.96"
    String spades_container = "microbiomedata/spades:3.15.0"
    Boolean paired = true

    call stage {
        input: container = bbtools_container,
               input_file = input_file
    }

    call bbcms {
        input: input_files= stage.assembly_input,
               container = bbtools_container,
               memory = memory,
               paired = paired
    }
    call assy {
        input: infile1 = bbcms.out1,
               infile2 = bbcms.out2,
               container = spades_container,
               threads = threads,
               paired = paired
    }
    call create_agp {
        input: scaffolds_in = assy.out,
               container = bbtools_container,
               rename_contig_prefix = rename_contig_prefix,
               memory = memory
    }
    call read_mapping_pairs {
        input: reads = stage.assembly_input,
               ref = create_agp.outcontigs,
               container = bbtools_container,
               memory = memory,
               threads = threads,
               paired = paired
    }

    call finish_asm {
        input: proj = proj,
               container = "microbiomedata/workflowmeta:1.1.1",
               input_file = stage.assembly_input,
               fasta = create_agp.outcontigs,
               scaffold = create_agp.outscaffolds,
               agp = create_agp.outagp,
               bam = read_mapping_pairs.outbamfile,
               samgz = read_mapping_pairs.outsamfile,
               covstats = read_mapping_pairs.outcovfile,
               asmstats = create_agp.outstats
    }
 
    output {
        File contig = finish_asm.outcontigs
        File scaffold = finish_asm.outscaffolds
        File agp = finish_asm.outagp
        File bam = finish_asm.outbam
        File samgz = finish_asm.outsamgz
        File covstats = finish_asm.outcovstats
        File asmstats = finish_asm.outasmstats
        Object stats = finish_asm.stats
        String start = stage.start
        String end = finish_asm.end
    }
 
    meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.0"
    }

}

task stage {
    String container
    String input_file
    String? memory = "4G"
    String target = "staged.fastq.gz"
    String output1 = "input.left.fastq.gz"
    String output2 = "input.right.fastq.gz"

    command <<<
        set -e
        if [ $( echo ${input_file}|egrep -c "https*:") -gt 0 ] ; then
            wget ${input_file} -O ${target}
        else
            ln ${input_file} ${target} || cp ${input_file} ${target}
        fi
 
         reformat.sh -Xmx${default="10G" memory} in=${target} out1=${output1} out2=${output2}    
        # Capture the start time
        date --iso-8601=seconds > start.txt

    >>>

    output{
        Array[File] assembly_input = [output1, output2]
        String start = read_string("start.txt")
    }
    runtime {
        cpu:  2
        maxRetries: 1
        docker: container
    }
}

task finish_asm {
    Array[File] input_file
    File fasta
    File scaffold
    File? agp
    File bam
    File? samgz
    File? covstats
    File asmstats
    String container
    String proj
    String prefix = sub(proj, ":", "_")
    String orig_prefix = "scaffold"
    String sed = "s/${orig_prefix}_/${proj}_/g"
 
    command<<<

        set -e
        date --iso-8601=seconds > end.txt

        ##RE-ID
        cat ${fasta} | sed ${sed} > ${prefix}_contigs.fna
        cat ${scaffold} | sed ${sed} > ${prefix}_scaffolds.fna
        cat ${covstats} | sed ${sed} > ${prefix}_covstats.txt
        cat ${agp} | sed ${sed} > ${prefix}_assembly.agp

        ## Bam file     
        samtools view -h ${bam} | sed ${sed} | \
           samtools view -hb -o ${prefix}_pairedMapped_sorted.bam
        ## Sam.gz file
        samtools view -h ${samgz} | sed ${sed} | \
           gzip -c - > ${prefix}_pairedMapped.sam.gz

        # Fix field names
        cat ${asmstats} |tr 'A-Z' 'a-z'| \
           jq 'del(.filename)' > stats.json

    >>>
    output {
        File outcontigs = "${prefix}_contigs.fna"
        File outscaffolds = "${prefix}_scaffolds.fna"
        File outagp = "${prefix}_assembly.agp"
        File outbam = "${prefix}_pairedMapped_sorted.bam"
        File outsamgz = "${prefix}_pairedMapped.sam.gz"
        File outcovstats = "${prefix}_covstats.txt"
        File outasmstats = "stats.json"
        Object stats = read_json("stats.json")
        String end = read_string("end.txt")
    }

    runtime {
        docker: container
        memory: "1 GiB"
        cpu:  1
    }
}

task read_mapping_pairs{
    Array[File] reads
    File ref
    String container
    String? memory
    String? threads
    Boolean paired = true
    String bbmap_interleaved_flag = if paired then 'interleaved=true' else 'interleaved=false'
    String filename_unsorted = "pairedMapped.bam"
    String filename_outsam = "pairedMapped.sam.gz"
    String filename_sorted = "pairedMapped_sorted.bam"
    String filename_sorted_idx = "pairedMapped_sorted.bam.bai"
    String filename_bamscript = "to_bam.sh"
    String filename_cov = "covstats.txt"
    String system_cpu = "$(grep \"model name\" /proc/cpuinfo | wc -l)"
    String jvm_threads = select_first([threads,system_cpu])

    runtime {
        docker: container
        memory: "120 GiB"
        cpu:  16
        maxRetries: 1
    }
    command{
        set -eo pipefail
        if [[ ${reads[0]}  == *.gz ]] ; then
             cat ${sep=" " reads} > infile.fastq.gz
             export mapping_input="infile.fastq.gz"
        fi
        if [[ ${reads[0]}  == *.fastq ]] ; then
             cat ${sep=" " reads} > infile.fastq
             export mapping_input="infile.fastq"
        fi
        bbmap.sh -Xmx${default="105G" memory} threads=${jvm_threads} nodisk=true ${bbmap_interleaved_flag} ambiguous=random in=$mapping_input ref=${ref} out=${filename_unsorted} covstats=${filename_cov} bamscript=${filename_bamscript}
        samtools sort -m100M -@ ${jvm_threads} ${filename_unsorted} -o ${filename_sorted}
        samtools index ${filename_sorted}
        reformat.sh -Xmx${default="105G" memory} in=${filename_unsorted} out=${filename_outsam} overwrite=true
        ln ${filename_cov} mapping_stats.txt
        rm $mapping_input
    }
    output{
        File outbamfile = filename_sorted
        File outbamfileidx = filename_sorted_idx
        File outcovfile = filename_cov
        File outsamfile = filename_outsam
    }
}

task create_agp {
    File scaffolds_in
    String? memory
    String container
    String rename_contig_prefix
    String prefix = "assembly"
    String filename_contigs = "${prefix}_contigs.fna"
    String filename_scaffolds = "${prefix}_scaffolds.fna"
    String filename_agp = "${prefix}.agp"
    String filename_legend = "${prefix}_scaffolds.legend"
    runtime {
        docker: container
        memory: "120 GiB"
        cpu:  16
    }
    command{
        fungalrelease.sh -Xmx${default="105G" memory} in=${scaffolds_in} out=${filename_scaffolds} outc=${filename_contigs} agp=${filename_agp} legend=${filename_legend} mincontig=200 minscaf=200 sortscaffolds=t sortcontigs=t overwrite=t
        if [ "${rename_contig_prefix}" != "scaffold" ]; then
            sed -i 's/scaffold/${rename_contig_prefix}_scf/g' ${filename_contigs} ${filename_scaffolds} ${filename_agp} ${filename_legend}
        fi
        bbstats.sh format=8 in=${filename_scaffolds} out=stats.json
        sed -i 's/l_gt50k/l_gt50K/g' stats.json
    }
    output{
        File outcontigs = filename_contigs
        File outscaffolds = filename_scaffolds
        File outagp = filename_agp
        File outstats = "stats.json"
        File outlegend = filename_legend
    }
}

task assy {
    File infile1
    File infile2
    String container
    String? threads
    String outprefix = "spades3"
    String filename_outfile = "${outprefix}/scaffolds.fasta"
    String filename_spadeslog  = "${outprefix}/spades.log"
    String system_cpu = "$(grep \"model name\" /proc/cpuinfo | wc -l)"
    String spades_cpu = select_first([threads,system_cpu])
    Boolean paired = true
    runtime {
        docker: container
        memory: "120 GiB"
        cpu:  16
    }
    command{
        set -eo pipefail
        if ${paired}; then
            spades.py -m 2000 -o ${outprefix} --only-assembler -k 33,55,77,99,127  --meta -t ${spades_cpu} -1 ${infile1} -2 ${infile2}
        else
            spades.py -m 2000 -o ${outprefix} --only-assembler -k 33,55,77,99,127 -t ${spades_cpu} -s ${infile1}
        fi
    }
    output {
        File out = filename_outfile
        File outlog = filename_spadeslog
    }
}

task bbcms {
    Array[File] input_files
    String container
    String? memory
    Boolean paired = true

    String filename_outfile = "input.corr.fastq.gz"
    String filename_outfile1 = "input.corr.left.fastq.gz"
    String filename_outfile2 = "input.corr.right.fastq.gz"
    String filename_readlen = "readlen.txt"
    String filename_outlog = "stdout.log"
    String filename_errlog = "stderr.log"
    String filename_kmerfile = "unique31mer.txt"
    String filename_counts = "counts.metadata.json"
    runtime {
        docker: container
        memory: "120 GiB"
        cpu:  16
    }

    command {
        set -eo pipefail
        if file --mime -b ${input_files[0]} | grep gzip; then
             cat ${sep=" " input_files} > infile.fastq.gz
             export bbcms_input="infile.fastq.gz"
        fi
        if file --mime -b ${input_files[0]} | grep plain; then
             cat ${sep=" " input_files} > infile.fastq
             export bbcms_input="infile.fastq"
        fi
        bbcms.sh -Xmx${default="105G" memory} metadatafile=${filename_counts} mincount=2 highcountfraction=0.6 in=$bbcms_input out=${filename_outfile} > >(tee -a ${filename_outlog}) 2> >(tee -a ${filename_errlog} >&2) && grep Unique ${filename_errlog} | rev |  cut -f 1 | rev  > ${filename_kmerfile}
        if ${paired}; then
            reformat.sh -Xmx${default="105G" memory} in=${filename_outfile} out1=${filename_outfile1} out2=${filename_outfile2}
        fi
        readlength.sh -Xmx${default="105G" memory} in=${filename_outfile} out=${filename_readlen}
        rm $bbcms_input
    }
    output {
        File out = filename_outfile
        File out1 = if paired then filename_outfile1 else filename_outfile
        File out2 = if paired then filename_outfile2 else filename_outfile
        File outreadlen = filename_readlen
        File stdout = filename_outlog
        File stderr = filename_errlog
        File outcounts = filename_counts
        File outkmer = filename_kmerfile
    }
}

