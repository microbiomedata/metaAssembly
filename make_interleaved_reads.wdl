workflow make_interleaved_reads {
	Array[File] input_files
	String output_file = "interleaved.fastq.gz"
	
	call interleave_reads {
		input: input_files = input_files, output_file = output_file
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

	Array[File] input_files
	String output_file = "interleaved.fastq.gz"
	
	command <<<
	    if file --mime -b ${input_files[0]} | grep gzip; then 
			paste <(gunzip -c ${input_files[0]} | paste - - - -) <(gunzip -c ${input_files[1]} | paste - - - -) | tr '\t' '\n' | gzip -c > ${output_file}
		else
			paste <(cat ${input_files[0]} | paste - - - -) <(cat ${input_files[1]} | paste - - - -) | tr '\t' '\n' | gzip -c > ${output_file}
		fi

	>>>
	
	runtime {
		#docker: "my_image"
	}
	
	output {
		File out_fastq = output_file
	}
}