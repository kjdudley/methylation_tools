def import_contig_list(text_file):
    text_file_handle = open(text_file, "r")
    contig_list = []
    for i in text_file_handle:
        contig_name = i.strip()
        contig_list.append(contig_name)

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

def import_intergene_coordinates(gene_annotation, contig_sequence_dict):

    intergene_annotation = {}
    for i in gene_annotation:
        start = gene_annotation[i][0]
        end = gene_annotation[i][1]
        intergene_start = []
        intergene_end = []
        if len(start) > 0:
            for j in range(len(start)):
                if j == 0:
                    if int(start[j]) > 1: # therefore must be a 5-prime intergenic region
                        print("start intergenic (5-prime):", 1, "end intergenic (5-prime):", int(start[j])-1)
                        intergene_start.append(1)
                        intergene_end.append(int(start[j])-1)
                    print ("start gene (first gene):", start[j], "end gene (first gene):", end[j])
                    if len(start) > 1: # therefore must be more than 1 gene on contig
                        if int(start[j+1]) > int(end[j])+1: # check that there is actually an intergenic region between genes
                            print("start intergenic:", int(end[j])+1, "end intergenic:", int(start[j+1])-1)
                            intergene_start.append(int(end[j])+1)
                            intergene_end.append(int(start[j+1])-1)
                    if len(start) == 1: # therefore must be only 1 gene on contig
                        if len(contig_sequence_dict[i]) > int(end[j]):
                            print("start intergenic (3-prime):", int(end[j])+1, "end intergenic (3-prime):", len(contig_sequence_dict[i]))
                            intergene_start.append(int(end[j])+1)
                            intergene_end.append(len(contig_sequence_dict[i]))
                elif j == len(start)-1:
                    print("start gene (last gene):", start[j], "end gene (last gene):", end[j])
                    if len(contig_sequence_dict[i]) > int(end[j]):
                        print("start intergenic (3-prime):", int(end[j])+1, "end intergenic (3-prime):", len(contig_sequence_dict[i]))
                        intergene_start.append(int(end[j])+1)
                        intergene_end.append(len(contig_sequence_dict[i]))
                elif len(start) > 1:
                    print("start gene (middle gene):", start[j], "end gene (middle gene):", end[j])
                    if len(contig_sequence_dict[i]) > int(end[j]):
                        if int(start[j+1]) > int(end[j])+1: # check that there is actually an intergenic region between genes
                            print("start intergenic (middle):", int(end[j])+1, "end intergenic (middle):", int(start[j+1])-1)
                            intergene_start.append(int(end[j])+1)
                            intergene_end.append(int(start[j+1])-1)
        else:
            print("start intergenic (no genes):", 1, "end intergenic (no genes):", len(contig_sequence_dict[i]))
            intergene_start.append(1)
            intergene_end.append(len(contig_sequence_dict[i]))
        intergene_annotation[i] = [intergene_start, intergene_end]
    return intergene_annotation

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

    contig_count = 0

    for i in contig_list:
        contig_count += 1
        print("Processing contig", contig_count)
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

def contig_to_gene(contig_cpg_coordinate_dict, gene_annotation):

    contig_gene_cpg_coordinate_methylation_dict = {}

    for i in contig_cpg_coordinate_dict:
        print(i)
        gene_cpg_coordinate_dict = {}
        for j in range(len(gene_annotation[i][0])):
            start = gene_annotation[i][0][j]
            end = gene_annotation[i][1][j]
            strand = gene_annotation[i][2][j]
            coordinate = []
            percent = []
            print("Gene start:", gene_annotation[i][0][j], "Gene end:", gene_annotation[i][1][j], "Gene strand:", gene_annotation[i][2][j])
            for k in contig_cpg_coordinate_dict[i]:
                if k >= int(gene_annotation[i][0][j]) and k <= int(gene_annotation[i][1][j]):
                    coordinate.append(k)
                    percent.append(contig_cpg_percent_methylation_dict[i][contig_cpg_coordinate_dict[i].index(k)])
                    print("coordinate:", k, "percent:", contig_cpg_percent_methylation_dict[i][contig_cpg_coordinate_dict[i].index(k)])
            gene_cpg_coordinate_dict[i+"."+start] = [int(start), int(end), strand, coordinate, percent]
        contig_gene_cpg_coordinate_methylation_dict[i] = gene_cpg_coordinate_dict
    return contig_gene_cpg_coordinate_methylation_dict 

def gene_coordinate_to_position(contig_gene_cpg_coordinate_methylation_dict):

    contig_gene_cpg_position_methylation_dict = {}

    for i in contig_gene_cpg_coordinate_methylation_dict:
        gene_dict = {}
        for j in contig_gene_cpg_coordinate_methylation_dict[i]:
            print(contig_gene_cpg_coordinate_methylation_dict[i][j][0], contig_gene_cpg_coordinate_methylation_dict[i][j][1], contig_gene_cpg_coordinate_methylation_dict[i][j][2])
            strand = contig_gene_cpg_coordinate_methylation_dict[i][j][2]
            start = contig_gene_cpg_coordinate_methylation_dict[i][j][0]
            end = contig_gene_cpg_coordinate_methylation_dict[i][j][1]
            length = end-start
            print(length)
            cpg_position = []
            cpg_percent = []
            for k in range(len(contig_gene_cpg_coordinate_methylation_dict[i][j][3])): # scroll through coordinates for given gene
                if strand == "+":
                    position = contig_gene_cpg_coordinate_methylation_dict[i][j][3][k]-start
                    percent = contig_gene_cpg_coordinate_methylation_dict[i][j][4][k]
                    cpg_position.append(position)
                    cpg_percent.append(percent)
                elif strand == "-":
                    #print("position =", length, "- (", contig_gene_cpg_coordinate_methylation_dict[i][j][3][k], "-", start, ")")
                    position = length-(contig_gene_cpg_coordinate_methylation_dict[i][j][3][k]-start) # flip positions around
                    percent = contig_gene_cpg_coordinate_methylation_dict[i][j][4][k]
                    cpg_position.append(position)
                    cpg_percent.append(percent)
                else:
                    print("Error!")
            gene_dict[j] = [strand, start, end, cpg_position, cpg_percent]
        contig_gene_cpg_position_methylation_dict[i] = gene_dict
    return contig_gene_cpg_position_methylation_dict
