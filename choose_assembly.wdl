version 1.0
import "jgi_assembly.wdl" as srma
import "make_interleaved_WDL/make_interleaved_reads.wdl" as int
import "jgi_meta_wdl/metagenome_improved/metaflye.wdl" as lrma

workflow jgi_metaAssembly{
    input {  
        Boolean shortRead
        String proj
        String prefix=sub(proj, ":", "_")
        # shortReads parameters
        String? memory
        String? threads
        # longReads parameters
        Array[File] input_fastq
        String flye_container = "staphb/flye:2.9.2"
        String flye_parameters = "--meta -o flye -t 32 --pacbio-hifi"
        String smrtlink_container = "bryce911/smrtlink:12.0.0.177059"
        String racon_container = "staphb/racon:1.4.20"
        String minimap2_container = "staphb/minimap2:2.25"
        String minimap2_parameters = "-a -x map-hifi -t 32"
        String samtools_container = "staphb/samtools:1.18"
        String bbtools_container = "microbiomedata/bbtools:38.96"
    }


    if (shortRead) {
    	if (length(input_fastq) > 1){
        	call int.make_interleaved_reads{
			input:
			input_files = input_fastq,
            container = bbtools_container
       		}
    	}
        call srma.jgi_metaASM{
            input:
            memory = memory,
            threads = threads,
            input_file = if length(input_fastq) > 1 then make_interleaved_reads.interleaved_fastq else input_fastq[0],
            proj = proj,
            bbtools_container = bbtools_container

        }
        
    }
    if (!shortRead) {
        call lrma.metaflye{
            input:
            proj = proj,
            input_fastq = input_fastq,
            flye_container = flye_container,
            flye_parameters = flye_parameters,
            smrtlink_container = smrtlink_container, 
            racon_container = racon_container,
            minimap2_container = minimap2_container,
            minimap2_parameters = minimap2_parameters, 
            samtools_container = samtools_container,
            bbtools_container = bbtools_container
        }
        call finish_lrasm{
            input: 
            proj = proj,
            prefix = prefix,
            container = bbtools_container,
            contigs = metaflye.final_contigs,
            bam = metaflye.final_bam, 
            scaffolds = metaflye.final_scaffolds,
            agp = metaflye.final_agp,
            legend = metaflye.final_legend,
            basecov = metaflye.final_basecov,
            sam = metaflye.final_sam,    
            output_file = metaflye.final_output_file,
            stats = metaflye.final_stats,
            summary_stats = metaflye.final_summary_stats,
            pileup_out = metaflye.final_pileup_out
        }
    }
    output {
        # long reads output
        File? lr_contigs = finish_lrasm.final_contigs
        File? lr_bam = finish_lrasm.final_bam 
        File? lr_scaffolds = finish_lrasm.final_scaffolds
        File? lr_agp = finish_lrasm.final_agp
        File? lr_legend = finish_lrasm.final_legend
        File? lr_basecov = finish_lrasm.final_basecov
        File? lr_sam = finish_lrasm.final_sam    
        File? lr_output_file = finish_lrasm.final_output_file
        File? lr_stats = finish_lrasm.final_stats
        File? lr_summary_stats = finish_lrasm.final_summary_stats
        File? lr_pileup_out = finish_lrasm.final_pileup_out
        File? lr_asminfo = metaflye.asminfo

        # short reads output
        File? sr_contig=jgi_metaASM.contig
        File? sr_scaffold=jgi_metaASM.scaffold
        File? sr_agp=jgi_metaASM.agp
        File? sr_bam=jgi_metaASM.bam
        File? sr_samgz=jgi_metaASM.samgz
        File? sr_covstats=jgi_metaASM.covstats
        File? sr_asmstats=jgi_metaASM.asmstats
        File? sr_asminfo=jgi_metaASM.asminfo

          
    }
}


task finish_lrasm {
    input{
    File contigs
    File bam
    File scaffolds
    File agp
    File legend
    File basecov
    File sam
    File output_file
    File stats
    File summary_stats
    File pileup_out
    String container
    String proj
    String prefix 
    String orig_prefix="scaffold"
    String sed="s/~{orig_prefix}_/~{proj}_/g"
    # String start
    }
    command<<<

        set -oeu pipefail
        end=`date --iso-8601=seconds`
        ln ~{output_file} ~{prefix}_read_count_report.txt
        ln ~{stats} ~{prefix}_contigs.sam.stats
        ln ~{summary_stats} ~{prefix}_summary.stats

        ##RE-ID
        cat ~{contigs} | sed ~{sed} > ~{prefix}_contigs.fna
        cat ~{scaffolds} | sed ~{sed} > ~{prefix}_scaffolds.fna
        cat ~{agp} | sed ~{sed} > ~{prefix}_assembly.agp
        cat ~{legend} | sed ~{sed} > ~{prefix}_assembly.legend
        cat ~{basecov} | sed ~{sed} > ~{prefix}_contigs.sorted.bam.pileup.basecov
        cat ~{pileup_out} | sed ~{sed} > ~{prefix}_contigs.sorted.bam.pileup.out

       ## Bam file     
       samtools view -h ~{bam} | sed ~{sed} | \
          samtools view -hb -o ~{prefix}_pairedMapped_sorted.bam
       ## Sam.gz file
       samtools view -h ~{sam} | sed ~{sed} | \
          gzip -c - > ~{prefix}_pairedMapped.sam.gz

    >>>
    output {
        File final_contigs = "~{prefix}_contigs.fna"
        File final_bam = "~{prefix}_pairedMapped_sorted.bam"
        File final_scaffolds = "~{prefix}_scaffolds.fna"
        File final_agp = "~{prefix}_assembly.agp"
        File final_legend = "~{prefix}_assembly.legend"
        File final_basecov = "~{prefix}_contigs.sorted.bam.pileup.basecov"
        File final_sam = "~{prefix}_pairedMapped.sam.gz"  
        File final_output_file = "~{prefix}_read_count_report.txt"
        File final_stats = "~{prefix}_contigs.sam.stats"
        File final_summary_stats = "~{prefix}_summary.stats"
        File final_pileup_out = "~{prefix}_contigs.sorted.bam.pileup.out"
    }

    runtime {
        docker: container
        memory: "1 GiB"
        cpu:  1
    }
}
