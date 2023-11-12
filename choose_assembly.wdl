version 1.0
import "jgi_assembly.wdl" as srma
import "jgi_meta_wdl/metagenome_improved/metaflye.wdl" as lrma
import "make_interleaved_WDL/make_interleaved_reads.wdl" as int

workflow jgi_assembly{
    input {  
        Boolean shortRead
        # shortReads parameters
        String? memory
        String? threads
        String proj
        # longReads parameters
        Array[File] input_fastq
        String flye_container
        String flye_parameters
        String smrtlink_container
        String racon_container
        String minimap2_container
        String minimap2_parameters
        String samtools_container
        String bbtools_container
    }

    if (length(input_fastq) > 1){
        call int.make_interleaved_reads{
            input:
            input_files = input_fastq
        }
    }

    if (shortRead) {
        call srma.jgi_metaASM{
            input:
            memory = memory,
            threads = threads,
            input_file = if length(input_fastq) > 1 then make_interleaved_reads.interleaved_fastq else input_fastq[0],
            proj = proj
        }
        
    }
    if (!shortRead) {
        call lrma.metaflye{
            input:
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
    }
    output {
        File? final_contigs = metaflye.final_contigs
        File? final_bam = metaflye.final_bam 
        File? final_scaffolds = metaflye.final_scaffolds
        File? final_agp = metaflye.final_agp
        File? final_legend = metaflye.final_legend
        File? final_basecov = metaflye.final_basecov
        File? final_sam = metaflye.final_sam    
        File? final_output_file = metaflye.final_output_file
        File? final_stats = metaflye.final_stats
        File? final_summary_stats = metaflye.final_summary_stats
        File? final_pileup_out = metaflye.final_pileup_out

        File? contig=jgi_metaASM.contig
        File? scaffold=jgi_metaASM.scaffold
        File? agp=jgi_metaASM.agp
        File? bam=jgi_metaASM.bam
        File? samgz=jgi_metaASM.samgz
        File? covstats=jgi_metaASM.covstats
        File? asmstats=jgi_metaASM.asmstats
        File? asminfo=jgi_metaASM.asminfo
    }
}
