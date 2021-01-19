import "metagenome_assy.wdl" as metagenome_assy
import "mapping.wdl" as mapping

workflow metagenome_assembly_and_alignment {
    Array[File] input_files

    call metagenome_assy.metagenome_assy as assy {
        input: input_files=input_files
    }
    call mapping.mapping {
       input: input_files=input_files, input_reference=assy.final_contigs
    }
}
