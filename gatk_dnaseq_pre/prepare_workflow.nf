#! /usr/bin/env nextflow

ref_dict = file( params.fasta_ref.replace(".fasta",".dict") )
outpath = params.outpath

process CreateSequenceGroupingTSV {

    input:
    file ref_dict

    output:
    file 'sequence_grouping.txt' into sequence_grouping
    file 'sequence_grouping_with_unmapped.txt' into sequence_grouping_with_unmapped

    container "job-definition://gatk-sequencegrouping"

    script:
    """
    #!/usr/bin/env python3
with open("${ref_dict}", "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

# We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in
# some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"

# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += "\\t" + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += "\\n" + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]

# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
with open("sequence_grouping.txt","w") as tsv_file:
  tsv_file.write(tsv_string)
  tsv_file.close()

tsv_string += '\\n' + "unmapped"

with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
  tsv_file_with_unmapped.write(tsv_string)
  tsv_file_with_unmapped.close()
    """
}

process SplitSequenceGroupingTSV {
    input:
    file input_grouping from sequence_grouping

    output:
    file "sequence_grouping_*.interval_list"
    file "sequence_group_list.txt"

    container "job-definition://gatk-sequencegrouping"
    publishDir "${outpath}/sequence_groups/", mode: 'copy', overwrite: false

    script:
    """
    #!/usr/bin/env python3
with open("${input_grouping}", "r") as handle:
    with open("sequence_group_list.txt","w") as index_file:
        count = 0
        for tsv_string in handle:
            with open("sequence_grouping_%d.interval_list" % count,"w") as tsv_file:
                for line_part in tsv_string.split('\\t'):
                    tsv_file.write(line_part + "\\n")
                tsv_file.close()
                index_file.write("${outpath}/sequence_groups/sequence_grouping_%d.interval_list\\n" % count)
                count += 1
        index_file.close()
    """
}

process SplitSequenceGroupingUnmappedTSV {
    input:
    file input_grouping from sequence_grouping_with_unmapped

    output:
    file "sequence_grouping_with_unmapped_*.interval_list"
    file "sequence_grouping_with_unmapped_list.txt"

    container "job-definition://gatk-sequencegrouping"
    publishDir "${outpath}/sequence_groups_unmapped/", mode: 'copy', overwrite: false

    script:
    """
    #!/usr/bin/env python3
with open("${input_grouping}", "r") as handle:
    with open("sequence_grouping_with_unmapped_list.txt","w") as index_file:
        count = 0
        for tsv_string in handle:
            with open("sequence_grouping_with_unmapped_%d.interval_list" % count,"w") as tsv_file:
                for line_part in tsv_string.split('\\t'):
                    tsv_file.write(line_part + "\\n")
                tsv_file.close()
                index_file.write("${outpath}/sequence_groups_unmapped/sequence_grouping_with_unmapped_%d.interval_list\\n" % count)
                count += 1
        index_file.close()
    """
}
