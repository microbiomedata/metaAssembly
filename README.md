# The Metagenome Assembly Pipeline

## Summary
This workflow is developed by Brian Foster at JGI and original from his [repo](https://gitlab.com/bfoster1/wf_templates/tree/master/templates). It take paired-end reads runs error correction by bbcms (BBTools). The clean reads are assembled by MetaSpades. After assembly, the reads are mapped back to contigs by bbmap (BBTools) for coverage information.

## Running Workflow in Cromwell
We provide three ways to run the workflow.  
1. `CromwellJtmShifter/`: The Cromwell run in head node send tasks to jtm-task-managers which will manages the tasks running on a computer node and using Shifter to run applications. 
2. `SlurmCromwellShifter/`: The submit script will request a node and launch the Cromwell.  The Cromwell manages the workflow by using Shifter to run applications. 
3. `CromwellSlurmShifter/`: The Cromwell run in head node and manages the workflow by submitting each step of workflow to compute node where applications were ran by Shifter.

Description of the files in each sud-directory:
 - `.wdl` file: the WDL file for workflow definition
 - `.json` file: the example input for the workflow
 - `.conf` file: the conf file for running Cromwell.
 - `.sh` file: the shell script for running the example workflow

## The Docker image and Dockerfile can be found here

[microbiomedata/bbtools:38.44](https://hub.docker.com/r/microbiomedata/bbtools)

[microbiomedata/spades:3.13.0](https://hub.docker.com/r/microbiomedata/spades)


## Input files

1. fastq (illumina paired-end interleaved fastq)
    
2. contig prefix for fasta header
    
3. output path

```
{
  "jgi_meta.input_file":["/global/cfs/projectdirs/m3408/ficus/11809.7.220839.TCCTGAG-ACTGCAT.fastq.gz"],
  "jgi_meta.rename_contig_prefix":"503125_160870",
  "jgi_meta.outdir":"/global/cfs/projectdirs/m3408/aim2/metagenome/assembly/ficus/503125_160870"
}
```

## Output files

Below is a part list of all output files. The main assembly contigs output is in final_assembly/assembly.contigs.fasta.

```
	├── bbcms
	│   ├── berkeleylab-jgi-meta-60ade422cd4e
	│   ├── counts.metadata.json
	│   ├── input.corr.fastq.gz
	│   ├── input.corr.left.fastq.gz
	│   ├── input.corr.right.fastq.gz
	│   ├── readlen.txt
	│   └── unique31mer.txt
	├── final_assembly
	│   ├── assembly.agp
	│   ├── assembly.contigs.fasta
	│   ├── assembly.scaffolds.fasta
	│   └── assembly.scaffolds.legend
	├── mapping
	│   ├── covstats.txt
	│   ├── pairedMapped.bam
	│   ├── pairedMapped.sam.gz
	│   ├── pairedMapped_sorted.bam
	│   └── pairedMapped_sorted.bam.bai
	└── spades3
		├── assembly_graph.fastg
		├── assembly_graph_with_scaffolds.gfa
		├── contigs.fasta
		├── contigs.paths
		├── scaffolds.fasta
		└── scaffolds.paths	
```

