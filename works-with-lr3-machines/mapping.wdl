workflow mapping {
    Array[File] input_files
    File input_reference

    String bbtools_container="bryce911/bbtools:38.86"

    if (length(input_files) == 1 ){
       call mappingtask as single_run {
           input: reads=input_files[0], reference=input_reference, container=bbtools_container
       }
    }
    if (length(input_files) > 1 ){
       scatter (input_file in input_files) {
           call mappingtask as multi_run {
               input: reads=input_file, reference=input_reference, container=bbtools_container
           }
       }
    }
    call finalize_bams{
        input: insing=single_run.outbamfile, inmult=multi_run.outbamfile, reference=input_reference, container=bbtools_container
    }
    output {
        File final_outbam = finalize_bams.outbam
        File final_outsam = finalize_bams.outsam
        File final_outbamidx = finalize_bams.outbamidx
        File final_outcov = finalize_bams.outcov
        File final_outflagstat = finalize_bams.outflagstat
    }
}

task finalize_bams {
    File? insing
    Array[File]? inmult
    File reference
    String container

    String single = if(defined(insing)) then "1" else "0"
    String java="-Xmx50g"
    String filename_outsam="pairedMapped.sam.gz"
    String filename_sorted="pairedMapped_sorted.bam"
    String filename_sorted_idx="pairedMapped_sorted.bam.bai"
    String filename_cov="pairedMapped_sorted.bam.cov"
    String filename_flagstat="pairedMapped_sorted.bam.flagstat"    
    String dollar="$"

    command{
        SECONDS=0
        if [ ${single} == "1" ]
        then
                cp ${insing} ${filename_sorted}      
        else
                samtools merge ${filename_sorted} ${sep=" " inmult}
        fi
        samtools index ${filename_sorted}
        reformat.sh threads=${dollar}(nproc) ${java} in=${filename_sorted} out=${filename_outsam} overwrite=true
        pileup.sh in=${filename_sorted} out=${filename_cov} ref=${reference}
        samtools flagstat ${filename_sorted} 1>| ${filename_flagstat} 2>| ${filename_flagstat}.e          

        hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
        printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    }

    runtime {
	  docker: container
	  poolname: "bfoster_ma_jeff"
      node: 1
      nwpn: 1
      mem: "200G"
      time: "03:00:00"
      shared: 0
      partition: "lr3"
      constraint: "lr3_c32,jgi_m256"
      account: "lr_jgicloud"
      qos: "condo_jgicloud"
    }

    output{
        File outsam = filename_outsam
        File outbam = filename_sorted
        File outbamidx = filename_sorted_idx
        File outcov = filename_cov
        File outflagstat = filename_flagstat    
    }
}

task mappingtask {
    File reads
    File reference

    String container 
    String filename_unsorted="pairedMapped.bam"
    String filename_sorted="pairedMapped_sorted.bam"
    String dollar="$"
    
    command{
      SECONDS=0
      bbmap.sh threads=${dollar}(nproc)  nodisk=true \
      interleaved=true ambiguous=random rgid=filename \
      in=${reads} ref=${reference} out=${filename_unsorted}

      samtools sort -m200M -@ ${dollar}(nproc) ${filename_unsorted} -o ${filename_sorted}

      hrs=$(( SECONDS/3600 )); mins=$(( (SECONDS-hrs*3600)/60)); secs=$(( SECONDS-hrs*3600-mins*60 ))
      printf 'Time spent: %02d:%02d:%02d\n' $hrs $mins $secs
    }

    runtime {
	  docker: container
	  poolname: "bfoster_ma_jeff"
      node: 1
      nwpn: 1
      mem: "200G"
      time: "03:00:00"
      shared: 0
      partition: "lr3"
      constraint: "lr3_c32,jgi_m256"
      account: "lr_jgicloud"
      qos: "condo_jgicloud"
    }

    output{
      File outbamfile = filename_sorted
   }
}
