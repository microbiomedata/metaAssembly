# Make interleaved workflows for QC, etc.
version 1.0
workflow make_interleaved_reads {
	input{
	Array[String] input_files
	String output_file = "interleaved.fastq.gz"
	String container="microbiomedata/bbtools:38.96"
	}
	call interleave_reads {
		input: 
		input_files = input_files, 
		output_file = output_file,
		container = container
	}
	output {
		File interleaved_fastq = interleave_reads.out_fastq
	}
	
	parameter_meta{
	    input_files: "two fastq files, [forwards.fastq, reverse.fastq]"
	    output_file: "output interleaved fastq file"
	}
	meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.0"
    }
}

task interleave_reads{
	input{
	Array[String] input_files
	String output_file = "interleaved.fastq.gz"
	String container
	}
	command <<<
		set -euo pipefail
		if file --mime -b ~{input_files[0]} | grep gzip; then 
			paste <(gunzip -c ~{input_files[0]} | paste - - - -) <(gunzip -c ~{input_files[1]} | paste - - - -) | tr '\t' '\n' | gzip -c > ~{output_file}
		else
			paste <(cat ~{input_files[0]} | paste - - - -) <(cat ~{input_files[1]} | paste - - - -) | tr '\t' '\n' | gzip -c > ~{output_file}
		fi

	>>>
	
	runtime {
		docker: container
		memory: "1 GiB"
        cpu:  1
	}
	
	output {
		File out_fastq = output_file
	}
}
