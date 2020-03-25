def import_contig_sequence(contig_list):

    '''a function that returns a dictionary of contig sequences given a name or list of names'''

    contig_input = input("Enter the path to the contig sequence file: ")
    contig_sequence_dict = {}
    contig_handle = open(contig_input, "r")
    contig_fasta = contig_handle.readlines()
    print(contig_fasta[0:3])
    contig_fasta.append(">\n") # add this so that last contig is included
    contig_handle.close()
    for i in contig_list:
        extract = 0
        contig_sequence = ""
        print(i)
        for j in contig_fasta:
            if extract == 1:
                if j[0] == ">":
                    extract = 0
                    contig_sequence_dict[i] = contig_sequence
                    print(contig_sequence[0:201])
                    break
                else:
                    line = j.strip()
                    contig_sequence = contig_sequence + line
            if i in j:
                print(j)
                extract = 1
    return contig_sequence_dict

def import_gff_data(contig_list, feature):

    gff_input = input("Enter the path to the contig sequence file: ")
    gff_handle = open(gff_input, "r")
    gff = gff_handle.readlines()
    feature_dict = {}

    for i in contig_list:
        feature_annotation = []
        print(i)
        start = []
        end = []
        strand = []
        for j in gff:
            if j[0] != "#":
                line = j.split("\t")
                #print(line[0])
                if line[0] == i:
                    if line[2] == feature:
                        start.append(line[3])
                        end.append(line[4])
                        strand.append(line[6])
        feature_annotation.append(start)
        feature_annotation.append(end)
        feature_annotation.append(strand)
        feature_dict[i] = feature_annotation
    return feature_dict

def get_coordinates(sequence_dict, subsequence):

    '''a function that loops through the nucleotides of a given sequence and returns the start and end coordinates'''

    coordinate_dict = {}
    for i in sequence_dict:
        start_coordinate_list = []
        end_coordinate_list = []
        sequence = sequence_dict[i]
        print(sequence[0:11])
        for j in range(len(sequence)):
            if j > len(sequence)-len(subsequence):
                break
            else:
                tracker = j
                start = j
                end = j+len(subsequence)
                foundsequence = ""
                for k in range(len(subsequence)):
                    if sequence[tracker] == subsequence[k]:
                        print("match:", sequence[tracker],"\t", subsequence[k])
                        foundsequence = foundsequence+sequence[tracker]
                        tracker += 1
                    else:
                        break
                #print("foundsequence:", foundsequence)
                if foundsequence == subsequence:
                    start_coordinate_list.append(start)
                    end_coordinate_list.append(end)
        coordinate_list = [start_coordinate_list, end_coordinate_list]
        coordinate_dict[i] = coordinate_list
    return coordinate_dict

def import_methylation_data(contig_list, sample_CX_data_file):

    '''import the context specific methylation data for a given sample and specified contigs'''

    methylated_coords_data_handle = open(sample_CX_data_file, "r")
    methylated_coords_data = methylated_coords_data_handle.readlines()
    methylated_coords_data_handle.close()
    methylated_coords_data_dict = {} # this is where the methylation data for a given sample will be stored (methylated and unmethylated counts)

    for i in contig_list:
        coordinate_list = []
        strand_list = []
        methylated_list = []
        unmethylated_list = []
        context_list = []
        trinucleotide_sequence_list = []
        for j in methylated_coords_data:
            split_data = j.split("\t")
            if split_data[0] == i:
                coordinate_list.append(int(split_data[1]))
                strand_list.append(split_data[2])
                methylated_list.append(int(split_data[3]))
                unmethylated_list.append(int(split_data[4]))
                context_list.append(split_data[5])
                trinucleotide_sequence_list.append(split_data[6].rstrip("\n"))
        data_list = [coordinate_list, strand_list, methylated_list, unmethylated_list, context_list, trinucleotide_sequence_list]
        methylated_coords_data_dict[i] = data_list

    return methylated_coords_data_dict
