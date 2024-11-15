Metagenome Assembly Workflow (v1.0.7)
=====================================

.. image:: workflow_assembly.png
   :scale: 60%
   :alt: Metagenome assembly workflow dependencies

Workflow Overview
-----------------

This workflow takes in paired-end Illumina short reads or PacBio long reads.

**Short Reads**:

In short reads, the workflow reformats the interleaved file into two FASTQ files for downstream tasks using bbcms (BBTools). The corrected reads are assembled using metaSPAdes. After assembly, the reads are mapped back to contigs by bbmap (BBTools) for coverage information. The `.wdl` (Workflow Description Language) file includes five tasks: *bbcms*, *assy*, *create_agp*, *read_mapping_pairs*, and *make_output*.

1. The *bbcms* task takes in interleaved FASTQ inputs, performs error correction, and reformats the interleaved FASTQ into two output FASTQ files for paired-end reads for the next tasks. 
2. The *assy* task performs metaSPAdes assembly.
3. Contigs and Scaffolds (output of metaSPAdes) are processed by the *create_agp* task to rename the FASTA header and generate an `AGP format <https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/>`_ which describes the assembly.
4. The *read_mapping_pairs* task maps reads back to the final assembly to generate coverage information.
5. The final *make_output* task collects all output files into the specified directory.

**Long Reads**:

In long reads, the workflow uses Flye for assembly, pbmm2 for alignment, Racon for polishing, and minimap2 for read mapping and coverage analysis. The :literal:`.wdl` (Workflow Description Language) file includes six tasks: *combine_fastq*, *assy*, *racon*, *format_assembly*, *map*, and *make_info_file*.

1. The *combine_fastq* task combines the input FASTQ files into a single FASTQ file, which is used as input for polishing and mapping tasks.
2. The *assy* task takes in the input FASTQ files and performs assembly using Flye.
3. The *racon* task cleans up the assembled contigs through two rounds of error correction using :literal:`pbmm2` and :literal:`Racon`.
4. The *format_assembly* task formats the polished assembly using BBTools' :literal:`fungalrelease.sh`, creating release-ready scaffolds and contigs, along with an `AGP format <https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/>`_ file and a legend file that describes the assembly.
5. The *map* task maps the input reads back to the final assembly using minimap2 to generate coverage data.
6. The final *make_info_file* task produces a summary file documenting tool versions, parameters, memory usage, and Docker containers used throughout the workflow.


Workflow Availability
---------------------

The workflow from GitHub uses all the listed Docker images to run all third-party tools.  

The workflow is available on GitHub: `https://github.com/microbiomedata/metaAssembly`  

The corresponding Docker images are available on DockerHub:

- `https://hub.docker.com/r/microbiomedata/spades`
- `https://hub.docker.com/r/microbiomedata/bbtools`
- `https://hub.docker.com/r/staphb/flye`
- `https://hub.docker.com/r/staphb/racon`
- `https://hub.docker.com/r/staphb/minimap2`
- `https://hub.docker.com/r/staphb/samtools`

Requirements for Execution
--------------------------

(Recommendations are in **bold**)  

- WDL-capable Workflow Execution Tool (**Cromwell**)
- Container Runtime that can load Docker images (**Docker v2.1.0.3 or higher**) 

Hardware Requirements
---------------------

**Memory: >40 GB RAM**

The memory requirement depends on the input complexity. Here is a simple estimation equation for the memory required based on kmers in the input file::

    predicted_mem = (kmers * 2.962e-08 + 1.630e+01) * 1.1 (in GB)

.. note::

   The kmers variable for the equation above can be obtained using the `kmercountmulti.sh` script from BBTools.

   Example command:

   ::

       kmercountmulti.sh -k=31 in=your.read.fq.gz

Workflow Dependencies
---------------------

Third-party software: (This is included in the Docker image.)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `metaSPAdes v4.0.0 <https://cab.spbu.ru/software/spades/>`_ (License: `GPLv2 <https://github.com/ablab/spades/blob/spades_3.15.0/assembler/GPLv2.txt>`_)
- `BBTools v39.03 <https://jgi.doe.gov/data-and-tools/bbtools/>`_ (License: `BSD-3-Clause-LBNL <https://bitbucket.org/berkeleylab/jgi-bbtools/src/master/license.txt>`_)

Sample dataset(s)
-----------------

**Short Reads:**

- Small dataset: `Ecoli 10x (287M) <https://portal.nersc.gov/cfs/m3408/test_data/metaAssembly_small_test_data.tgz>`_ (Input/output included in tar.gz file)
- Large dataset: `Zymobiomics mock-community DNA control (22G) <https://portal.nersc.gov/cfs/m3408/test_data/metaAssembly_large_test_data.tgz>`_ (Input/output included in tar.gz file)
- Zymobiomics mock-community DNA control (`SRR7877884 <https://www.ebi.ac.uk/ena/browser/view/SRR7877884>`_). The original dataset is ~4 GB. For testing, a 10% subsample of the dataset is used: (`SRR7877884-int-0.1.fastq.gz <https://portal.nersc.gov/cfs/m3408/test_data/SRR7877884-int-0.1.fastq.gz>`_). This dataset is already interleaved. 

**Long Reads:**

Zymobiomics synthetic metagenome (`SRR13128014 <https://portal.nersc.gov/project/m3408//test_data/SRR13128014.pacbio.subsample.ccs.fastq.gz>`_) For testing we have subsampled the dataset, the original dataset is ~18GB.


Input
-----

A `JSON file <https://github.com/microbiomedata/metaAssembly/blob/master/input.json>`_ containing the following information:

1. The path to the input FASTQ file (Illumina paired-end interleaved FASTQ or PacBio paired-end interleaved FASTQ) (recommended: output of the Reads QC workflow).
2. Project name: nmdc:XXXXXX
3. Memory (optional) e.g., `"jgi_metaAssembly.memory": "105G"`
4. Threads (optional) e.g., `"jgi_metaAssembly.threads": "16"`
5. Whether the input is short reads (boolean)

Example input JSON for short reads::

    {
        "jgi_metaAssembly.input_files": ["https://portal.nersc.gov/project/m3408/test_data/smalltest.int.fastq.gz"],
        "jgi_metaAssembly.proj": "nmdc:XXXXXX",
        "jgi_metaAssembly.memory": "105G",
        "jgi_metaAssembly.threads": "16",
        "jgi_metaAssembly.shortRead": true
    }

Example input JSON for long reads::

    {
        "jgi_metaAssembly.input_files": ["/global/cfs/cdirs/m3408/www/test_data/SRR13128014.pacbio.subsample.ccs.fastq.gz"],
        "jgi_metaAssembly.proj": "nmdc:XXXXXX",
        "jgi_metaAssembly.memory": "105G",
        "jgi_metaAssembly.threads": "16",
        "jgi_metaAssembly.shortRead": false
    }

Output
------

The output directory will contain the following files for short reads::

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

The output directory will contain the following files for long reads::

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

Example output stats JSON file::

    {
       "scaffolds": 58,
       "contigs": 58,
       "scaf_bp": 28406,
       "contig_bp": 28406,
       "gap_pct": 0.00000,
       "scaf_N50": 21,
       "scaf_L50": 536,
       "ctg_N50": 21,
       "ctg_L50": 536,
       "scaf_N90": 49,
       "scaf_L90": 317,
       "ctg_N90": 49,
       "ctg_L90": 317,
       "scaf_logsum": 22.158,
       "scaf_powsum": 2.245,
       "ctg_logsum": 22.158,
       "ctg_powsum": 2.245,
       "asm_score": 0.000,
       "scaf_max": 1117,
       "ctg_max": 1117,
       "scaf_n_gt50K": 0,
       "scaf_l_gt50K": 0,
       "scaf_pct_gt50K": 0.0,
       "gc_avg": 0.39129,
       "gc_std": 0.03033
    }


The table provides all of the output directories, files, and their descriptions.


=================================================== ===================================================== ===============================================================
Directory                                           File Name                                             Description
=================================================== ===================================================== ===============================================================
**Short Reads**                                                                                           Short reads assembly output directory
/make_info_file                                     nmdc_XXXXXX_metaAsm.info                              Summary information about the short reads assembly process
/finish_asm                                         nmdc_XXXXXX_covstats.txt                              Coverage statistics for assembled contigs
/finish_asm                                         nmdc_XXXXXX_contigs.fna                               Final contig sequences in FASTA format
/finish_asm                                         nmdc_XXXXXX_bbcms.fastq.gz                            Error-corrected FASTQ file from bbcms
/finish_asm                                         nmdc_XXXXXX_scaffolds.fna                             Final scaffold sequences in FASTA format
/finish_asm                                         nmdc_XXXXXX_assembly.agp                              Assembly information in AGP format
/finish_asm                                         stats.json                                            Assembly statistics in JSON format
/finish_asm                                         nmdc_XXXXXX_pairedMapped.sam.gz                       SAM file with reads mapped back to assembly
/finish_asm                                         nmdc_XXXXXX_pairedMapped_sorted.bam                   Sorted BAM file with reads mapped back to assembly

**Long Reads**                                                                                            Long reads assembly output directory
/finish_lrasm                                        nmdc_XXXXXX_assembly.legend                          Mapping file from contig to scaffold names
/finish_lrasm                                        nmdc_XXXXXX_contigs.fna                              Final contig sequences in FASTA format
/finish_lrasm                                        nmdc_XXXXXX_pairedMapped_sorted.bam                  Sorted BAM file with reads mapped back to assembly
/finish_lrasm                                        nmdc_XXXXXX_read_count_report.txt                    Read count report for validation
/make_info_file                                      nmdc_XXXXXX_metaAsm.info                             Summary information about the long reads assembly process
/finish_lrasm                                        nmdc_XXXXXX_summary.stats                            Summary statistics for assembly
/finish_lrasm                                        nmdc_XXXXXX_scaffolds.fna                            Final scaffold sequences in FASTA format
/finish_lrasm                                        nmdc_XXXXXX_pairedMapped.sam.gz                      SAM file with reads mapped back to assembly
/finish_lrasm                                        stats.json                                           Assembly statistics in JSON format
/finish_lrasm                                        nmdc_XXXXXX_contigs.sam.stats                        SAM file statistics for contigs
/finish_lrasm                                        nmdc_XXXXXX_contigs.sorted.bam.pileup.basecov        Base coverage information for contigs
/finish_lrasm                                        nmdc_XXXXXX_assembly.agp                             Assembly information in AGP format
/finish_lrasm                                        nmdc_XXXXXX_contigs.sorted.bam.pileup.out            BAM file pileup output for contigs
=================================================== ===================================================== ===============================================================


Version History
---------------

- 1.0.7 (release date **11/14/24**; previous versions: 1.0.6)

Point of contact
----------------

- Original author: Brian Foster <bfoster@lbl.gov>

- Package maintainer: Chienchi Lo <chienchi@lanl.gov>
