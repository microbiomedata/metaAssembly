# The Metagenome Assembly Pipeline

## Summary
This workflow is developed by Brian Foster at JGI and original from his [repo](https://gitlab.com/bfoster1/wf_templates/tree/master/templates). It take paired-end reads runs error correction by bbcms (BBTools). The clean reads are assembled by MetaSpades. After assembly, the reads are mapped back to contigs by bbmap (BBTools) for coverage information.

## Running Workflow in Cromwell

Description of the files:
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

4. memory (optional) ex: "jgi_metaASM.memory": "105G"

5. threads (optional) ex: "jgi_metaASM.threads": "16"

```
{
  "jgi_metaASM.input_file":["/global/cfs/projectdirs/m3408/ficus/11809.7.220839.TCCTGAG-ACTGCAT.fastq.gz"],
  "jgi_metaASM.rename_contig_prefix":"503125_160870",
  "jgi_metaASM.outdir":"/global/cfs/projectdirs/m3408/aim2/metagenome/assembly/ficus/503125_160870",
  "jgi_metaASM.memory": "105G",
  "jgi_metaASM.threads": "16"
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
	│   ├── assembly_contigs.fasta
	│   ├── assembly_scaffolds.fasta
	│   └── assembly_scaffolds.legend
	├── mapping
	│   ├── covstats.txt (mapping_stats.txt)
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

