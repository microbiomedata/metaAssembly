workflow metagenome_assy {
    Array[File] input_files
    Boolean nersc = false

    String bbtools_container="bryce911/bbtools:38.86"
    String spades_container="bryce911/spades:3.14.1"

    call bbcms {
    	 input: reads_files=input_files, container=bbtools_container
    }
    call assy {
    	 input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container    
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    output {
        File final_contigs = create_agp.outcontigs
        File final_scaffolds = create_agp.outscaffolds
        File final_spades_log = assy.outlog
	    File final_readlen = bbcms.outreadlen
	    File final_counts = bbcms.outcounts
    }
}

task bbcms{
    Array[File] reads_files

    String container
    String single = if (length(reads_files) == 1 ) then "1" else "0"

    String bbcms_input = "bbcms.input.fastq.gz"
    String filename_counts="counts.metadata.json"
  
    String filename_outfile="input.corr.fastq.gz"
    String filename_outfile1="input.corr.left.fastq.gz"
    String filename_outfile2="input.corr.right.fastq.gz"
    
    String filename_readlen="readlen.txt"
    String filename_outlog="stdout.log"
    String filename_errlog="stderr.log"
    String filename_kmerfile="unique31mer.txt"

    String java="-Xmx20g"
    String dollar="$"

    command {
		SECONDS=0
		
        if [ ${single} == 0 ]
	    then
	        cat ${sep = " " reads_files } > ${bbcms_input}
	    else
	        ln -s ${reads_files[0]} ./${bbcms_input}
	    fi

		touch ${filename_readlen}
		readlength.sh -Xmx1g in=${bbcms_input} out=${filename_readlen} overwrite 

		bbcms.sh ${java} metadatafile=${filename_counts} mincount=2 highcountfraction=0.6 \
	    in=${bbcms_input} out1=${filename_outfile1} out2=${filename_outfile2} \
	    1> ${filename_outlog} 2> ${filename_errlog}

		reformat.sh  in1=${filename_outfile1} in2=${filename_outfile2}  \
	    out=${filename_outfile} 1>> ${filename_outlog} 2>> ${filename_errlog}

		hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
		printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs

     }

     runtime {
		docker: container
		time: "01:00:00"
		mem: "115G"
		poolname: "bfoster_ma_wdl"
		node: 1
		nwpn: 1
		shared: 0
     }

     output {
        File out = filename_outfile
        File out1 = filename_outfile1
        File out2 = filename_outfile2
        File outreadlen = filename_readlen
        File stdout = filename_outlog
        File stderr = filename_errlog
        File outcounts = filename_counts
     }

}

task assy{
    File infile1
    File infile2    

    String container

    String outprefix="spades3"
    String filename_outfile="${outprefix}/scaffolds.fasta"
    String filename_spadeslog ="${outprefix}/spades.log"
    String dollar="$"

    command{
	   SECONDS=0
       spades.py -m 2000 --tmp-dir ${dollar}PWD -o ${outprefix} --only-assembler -k 33,55,77,99,127 --meta -t ${dollar}(nproc) -1 ${infile1} -2 ${infile2}
	   hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
	   printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    }

    runtime {
	  docker: container
	  time: "01:00:00"
	  mem: "115G"
	  poolname: "bfoster_ma_wdl"
	  node: 1
	  nwpn: 1
	  shared: 0
    }

    output {
           File out = filename_outfile
           File outlog = filename_spadeslog
    }
}

task create_agp {
    File scaffolds_in
    String container
    String prefix="assembly"

    String filename_contigs="${prefix}.contigs.fasta"
    String filename_scaffolds="${prefix}.scaffolds.fasta"
    String filename_agp="${prefix}.agp"
    String filename_legend="${prefix}.scaffolds.legend"

    command{
	  SECONDS=0
      fungalrelease.sh -Xmx105g in=${scaffolds_in} out=${filename_scaffolds} \
      outc=${filename_contigs} agp=${filename_agp} legend=${filename_legend} \
      mincontig=200 minscaf=200 sortscaffolds=t sortcontigs=t overwrite=t

	  hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
	  printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    }

    runtime {
	  docker: container
	  time: "01:00:00"
	  mem: "115G"
	  poolname: "bfoster_ma_wdl"
	  node: 1
	  nwpn: 1
	  shared: 0
    }

    output{
        File outcontigs = filename_contigs
        File outscaffolds = filename_scaffolds
        File outagp = filename_agp
        File outlegend = filename_legend
    }
}
