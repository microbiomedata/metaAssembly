# The Metagenome Assembly Pipeline

## Summary
This workflow is developed by Brian Foster at JGI and original from his [repo](https://gitlab.com/bfoster1/wf_templates/tree/master/templates). It takes in paired-end Illumina short reads or PacBio long reads. 

In short reads, the workflow reformats the interleaved file into two FASTQ files for downstream tasks using bbcms (BBTools). The corrected reads are assembled using metaSPAdes. After assembly, the reads are mapped back to contigs by bbmap (BBTools) for coverage information. The `.wdl` (Workflow Description Language) file includes five tasks: *bbcms*, *assy*, *create_agp*, *read_mapping_pairs*, and *make_output*.

In long reads, the workflow uses Flye for assembly, pbmm2 for alignment, Racon for polishing, and minimap2 for read mapping and coverage analysis. The `.wdl` (Workflow Description Language) file includes six tasks: *combine_fastq*, *assy*, *racon*, *format_assembly*, *map*, and *make_info_file*.


## The Docker image and Dockerfile can be found here

[microbiomedata/bbtools:39.03](https://hub.docker.com/r/microbiomedata/bbtools)

[microbiomedata/spades:4.0.0](https://hub.docker.com/r/microbiomedata/spades)


## Input files

1. The path to the input FASTQ file (Illumina paired-end interleaved FASTQ or PacBio paired-end interleaved FASTQ) (recommended: output of the Reads QC workflow).
    
2. Project name, e.g. `nmdc:XXXXXX`
    
3. Memory (optional) e.g., `"jgi_metaAssembly.memory": "105G"`

4. Threads (optional) e.g., `"jgi_metaAssembly.threads": "16"`

5. Whether the input is short reads (boolean) 


```
{
        "jgi_metaAssembly.input_files": ["https://portal.nersc.gov/project/m3408/test_data/smalltest.int.fastq.gz"],
        "jgi_metaAssembly.proj": "nmdc:XXXXXX",
        "jgi_metaAssembly.memory": "105G",
        "jgi_metaAssembly.threads": "16",
        "jgi_metaAssembly.shortRead": true
}
```

## Output files

Below is a part list of all output files. The main assembly contigs output is in final_assembly/assembly.contigs.fasta.

```
# Short Reads
    output/
    ├── nmdc_XXXXXX_metaAsm.info
    ├── nmdc_XXXXXX_covstats.txt
    ├── nmdc_XXXXXX_contigs.fna
    ├── nmdc_XXXXXX_bbcms.fastq.gz
    ├── nmdc_XXXXXX_scaffolds.fna
    ├── nmdc_XXXXXX_assembly.agp
    ├── stats.json
    ├── nmdc_XXXXXX_pairedMapped.sam.gz
    └── nmdc_XXXXXX_pairedMapped_sorted.bam
# Long Reads
    output/
    ├── nmdc_XXXXXX_assembly.legend
    ├── nmdc_XXXXXX_contigs.fna
    ├── nmdc_XXXXXX_pairedMapped_sorted.bam
    ├── nmdc_XXXXXX_read_count_report.txt
    ├── nmdc_XXXXXX_metaAsm.info
    ├── nmdc_XXXXXX_summary.stats
    ├── nmdc_XXXXXX_scaffolds.fna
    ├── nmdc_XXXXXX_pairedMapped.sam.gz
    ├── stats.json
    ├── nmdc_XXXXXX_contigs.sam.stats
    ├── nmdc_XXXXXX_contigs.sorted.bam.pileup.basecov
    ├── nmdc_XXXXXX_assembly.agp
    └── nmdc_XXXXXX_contigs.sorted.bam.pileup.out
```
## Link to Doc Site
Please refer [here](https://docs.microbiomedata.org/workflows/chapters/4_Metagenome_Assembly/) for more information.

