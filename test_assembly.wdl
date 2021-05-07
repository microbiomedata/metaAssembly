import "jgi_assembly.wdl" as jgi

workflow test_assembly {
  String bbtools_container="microbiomedata/bbtools:38.90"
  String spades_container="microbiomedata/spades:3.15.0"
  String validate_container="mbabinski17/comparejson:0.1"
  String rename_contig_prefix="scaffold"
  Float  uniquekmer=1000
  String? memory="60G"
  String? threads="8"
  String? outdir="/vol_b/nmdc_workflows/test_nmdc/metaAssembly/test_output"
  String  url="https://portal.nersc.gov/cfs/m3408/test_data/Ecoli_10x-int.fastq.gz"
  String  ref_json="https://raw.githubusercontent.com/microbiomedata/metaAssembly/master/test_output/small_test_stats.json"

  call prepare {
    input: container=bbtools_container,
	   url=url,
           ref_json=ref_json
  }
  call jgi.jgi_metaASM as asm {
    input: input_file=[prepare.fastq],
           bbtools_container=bbtools_container,
	   spades_container=spades_container,
           memory=memory,
           threads=threads,
	   uniquekmer=uniquekmer,
	   rename_contig_prefix=rename_contig_prefix,
	   outdir=outdir
  }
  call validate {
    input: container=validate_container,
           refjson=prepare.refjson,
           user_json=asm.final_asmstat
  }
}
task prepare {
   String container
   String ref_json
   String url
   command{
       wget -O "input.fastq.gz" ${url}
       wget -O "ref_json.json" ${ref_json}
   }
   output{
      File fastq = "input.fastq.gz"
      File refjson = "ref_json.json"
   }
   runtime {
     memory: "1 GiB"
     cpu:  2
     maxRetries: 1
     docker: container
   }
}
task validate {
   String container
   File refjson
   File user_json
   command {
       compare_json.py -i ${refjson} -f ${user_json}
   }
   output {
       Array[String] result = read_lines(stdout())
   }
   runtime {
     memory: "1 GiB"
     cpu:  1
     maxRetries: 1
     docker: container
   }
}
