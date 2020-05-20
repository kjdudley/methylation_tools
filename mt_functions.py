import time
import statistics

def import_contig_list(text_file_path):
    contig_list = []
    with open(text_file_path, "r") as text_file:
        for i in text_file:
            contig_name = i.strip()
            contig_list.append(contig_name)
    return contig_list

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

    """a function that returns a dictionary containg a list of start coordinates
    and a list of end coordinates for non genic regions (i.e. intergenic regions)"""

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

    '''a function that loops through the nucleotides of a given sequence and returns the start
    and end coordinates of any occurence of a subsequence'''

    coordinate_dict = {}
    for i in sequence_dict:
        start_coordinate_list = []
        end_coordinate_list = []
        sequence = sequence_dict[i]
        #print(sequence[0:11])
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

def import_methylation_data(contig_list, sample_CX_data_file, required_coverage):

    '''import the context specific methylation data for a given sample and specified contigs'''

    methylated_coords_data_handle = open(sample_CX_data_file, "r")
    methylated_coords_data = methylated_coords_data_handle.readlines()
    methylated_coords_data_handle.close()
    methylated_coords_data_dict = {} # this is where the methylation data for a given sample will be stored (methylated and unmethylated counts)

    contig_count = 0

    for i in contig_list:
        contig_count += 1
        #print("Processing contig", contig_count)
        coordinate_list = []
        strand_list = []
        methylated_list = []
        unmethylated_list = []
        context_list = []
        trinucleotide_sequence_list = []
        for j in methylated_coords_data:
            split_data = j.split("\t")
            if split_data[0] == i:
                if int(split_data[3]) + int(split_data[4]) >= required_coverage:
                    coordinate_list.append(int(split_data[1]))
                    strand_list.append(split_data[2])
                    methylated_list.append(int(split_data[3]))
                    unmethylated_list.append(int(split_data[4]))
                    context_list.append(split_data[5])
                    trinucleotide_sequence_list.append(split_data[6].rstrip("\n"))
        data_list = [coordinate_list, strand_list, methylated_list, unmethylated_list, context_list, trinucleotide_sequence_list]
        methylated_coords_data_dict[i] = data_list

    return methylated_coords_data_dict

def import_methylation_data2(contig_list, sample_CX_data_file, required_coverage):
    '''import the context specific methylation data for a given sample and specified contigs'''
    methylated_coords_data_handle = open(sample_CX_data_file, "r")
    methylated_coords_data = methylated_coords_data_handle.readlines()
    methylated_coords_data_handle.close()
    methylated_coords_data_dict = {}
    for i in contig_list:
        methylated_coords_data_dict[i] = [[],[],[],[],[],[]]
    icount = 0
    start = time.time()
    #print("starting data import:", start)
    #print("number of lines to import:", len(methylated_coords_data))
    hundreth = len(methylated_coords_data)//100
    increment = int(hundreth)
    percent = 0
    #print("***", percent, "percent complete ***")
    for i in methylated_coords_data:
        # print(icount)
        icount += 1
        split_data = i.split("\t")
        contig = split_data[0]
        if int(split_data[3]) + int(split_data[4]) >= required_coverage:
            methylated_coords_data_dict[contig][0] = methylated_coords_data_dict[contig][0] + [int(split_data[1])] # coordinate
            methylated_coords_data_dict[contig][1] = methylated_coords_data_dict[contig][1] + [(split_data[2])] # strand
            methylated_coords_data_dict[contig][2] = methylated_coords_data_dict[contig][2] + [(int(split_data[3]))] # methylated count
            methylated_coords_data_dict[contig][3] = methylated_coords_data_dict[contig][3] + [(int(split_data[4]))] # unmethylated count
            methylated_coords_data_dict[contig][4] = methylated_coords_data_dict[contig][4] + [(split_data[5])] # context
            methylated_coords_data_dict[contig][5] = methylated_coords_data_dict[contig][5] + [(split_data[6].rstrip("\n"))] # trinucleotide sequence
        if icount == increment:
            percent += int(100/100)
            now = time.time()
            #print("***", percent, "percent complete [", now-start, "seconds elapsed ] ***")
            increment += hundreth
    return methylated_coords_data_dict

def import_methylation_data3(input_file_path):
    '''import the context specific methylation data for a given sample and specified contigs'''
    contig_list = import_contig_list("/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest/contigs.txt")
    required_coverage = 20
    methylated_coords_data_handle = open(input_file_path, "r")
    methylated_coords_data = methylated_coords_data_handle.readlines()
    methylated_coords_data_handle.close()
    methylated_coords_data_dict = {} # this is where the methylation data for a given sample will be stored (methylated and unmethylated counts)
    log_file = open("log_file.txt", "a+")
    log_file.write('{0}/{1}/{2}\t{3}:{4}'.format(time.localtime()[2], time.localtime()[1], time.localtime()[0], time.localtime()[3], time.localtime()[4])+"\t*** Started processing "+ input_file_path+" ***\n")
    log_file.close()
    contig_count = 0
    for i in contig_list:
        contig_count += 1
        print("Processing contig", contig_count, "("+input_file_path+")")
        coordinate_list = []
        strand_list = []
        methylated_list = []
        unmethylated_list = []
        context_list = []
        trinucleotide_sequence_list = []
        for j in methylated_coords_data:
            split_data = j.split("\t")
            if split_data[0] == i:
                if int(split_data[3]) + int(split_data[4]) >= required_coverage: # required coverage
                    coordinate_list.append(int(split_data[1]))
                    strand_list.append(split_data[2])
                    methylated_list.append(int(split_data[3]))
                    unmethylated_list.append(int(split_data[4]))
                    context_list.append(split_data[5])
                    trinucleotide_sequence_list.append(split_data[6].rstrip("\n"))
        data_list = [coordinate_list, strand_list, methylated_list, unmethylated_list, context_list, trinucleotide_sequence_list]
        methylated_coords_data_dict[i] = data_list
    log_file = open("log_file.txt", "a+")
    log_file.write('{0}/{1}/{2}\t{3}:{4}'.format(time.localtime()[2], time.localtime()[1], time.localtime()[0], time.localtime()[3], time.localtime()[4])+"\t*** Finished processing "+ input_file_path+" ***\n")
    log_file.close()
    return methylated_coords_data_dict

def get_feature_stats(contig_cpg_coordinate_dict, contig_feature_coordinate_dict, contig_sequence_dict):
    """extract gene stats (i.e. size, cpg count, cpg percent, cpg coordinates)"""
    contig_feature_stats_dict = {}
    for i in contig_cpg_coordinate_dict:
        #print(i)
        number_of_cpg = []
        cpg_coord = []
        feature_cpg_coord_dict = {}
        size_of_feature = []
        feature_percent_cpg = [] # i.e proportion of feature made up of CpG dinucleotides
        for j in range(len(contig_feature_coordinate_dict[i][0])): # loop through features
            feature_start_coordinate = int(contig_feature_coordinate_dict[i][0][j]) # start coordinate
            feature_end_coordinate = int(contig_feature_coordinate_dict[i][1][j]) # end coordinate
            feature_cpg_count = 0
            feature_cpg_coord = []
            for k in range(len(contig_cpg_coordinate_dict[i][0])):
                if int(contig_cpg_coordinate_dict[i][0][k]) <= feature_end_coordinate:
                    if int(contig_cpg_coordinate_dict[i][0][k]) >= feature_start_coordinate:
                        feature_cpg_count += 1
                        feature_cpg_coord.append(contig_cpg_coordinate_dict[i][0][k])
                        #print(i,"\t",contig_cpg_coordinate_dict[i][0][k],"\t", contig_sequence_dict[i][contig_cpg_coordinate_dict[i][0][k]:contig_cpg_coordinate_dict[i][1][k]])
            number_of_cpg.append(feature_cpg_count)
            cpg_coord.append(feature_cpg_coord)
            feature_cpg_coord_dict[feature_start_coordinate] = feature_cpg_coord
            size_of_feature.append(feature_end_coordinate-(feature_start_coordinate-1))
            if feature_end_coordinate-feature_start_coordinate > 0:
                feature_percent_cpg.append(100*(feature_cpg_count/(feature_end_coordinate-feature_start_coordinate)))
            #print("feature_cpg_count:", feature_cpg_count)
        feature_stats = [feature_cpg_coord_dict, number_of_cpg, size_of_feature, feature_percent_cpg]
        contig_feature_stats_dict[i] = feature_stats
    return contig_feature_stats_dict

def get_sample_methylation_stats(sample_methylation_dict, required_coverage, context):
    contig_sample_methylation_stats_dict = {}
    contig_count = 0
    for i in sample_methylation_dict: # contig name
        contig_count += 1
        #print("Processing contig", contig_count, "Contig name:", i)
        contig_coordinate = []
        contig_percent_data = []
        contig_methylated_count_data = []
        contig_unmethylated_count_data = []
        contig_cumulative_methylated = 0
        contig_cumulative_unmethylated = 0
        contig_pass_coverage_threshold = 0
        contig_fail_coverage_threshold = 0
        for j in range(len(sample_methylation_dict[i][0])):
            if sample_methylation_dict[i][4][j] == context: # remove "and sample_methylation_dict[i][1][j] == '+'"
                coordinate = sample_methylation_dict[i][0][j]
                methylated = sample_methylation_dict[i][2][j]
                unmethylated = sample_methylation_dict[i][3][j]
                total = methylated+unmethylated
                context = sample_methylation_dict[i][4][j]
                if total >= required_coverage:
                    contig_coordinate.append(coordinate)
                    contig_methylated_count_data.append(methylated)
                    contig_unmethylated_count_data.append(unmethylated)
                    percent = 100*(methylated/total)
                    contig_percent_data.append(percent)
                    contig_cumulative_methylated += methylated
                    contig_cumulative_unmethylated += unmethylated
                    contig_pass_coverage_threshold += 1
                else:
                    contig_fail_coverage_threshold += 1
        contig_cpg_coverage_threshold_list = [contig_pass_coverage_threshold, contig_fail_coverage_threshold]
        contig_cpg_coordinate_list = contig_coordinate
        contig_cpg_methylated_count_list = contig_methylated_count_data
        contig_cpg_unmethylated_count_list = contig_unmethylated_count_data
        contig_cpg_percent_methylation_list = contig_percent_data
        contig_totals_data = [contig_cumulative_methylated, contig_cumulative_unmethylated]
        contig_cpg_counts_methylation_list = contig_totals_data
        sample_methylation_stats = [contig_cpg_coverage_threshold_list, contig_cpg_counts_methylation_list,  contig_cpg_coordinate_list, contig_cpg_percent_methylation_list, contig_cpg_methylated_count_list, contig_cpg_unmethylated_count_list]
        contig_sample_methylation_stats_dict[i] = sample_methylation_stats
#        if contig_cumulative_methylated+contig_cumulative_unmethylated > 0: # discounts contigs without any cytosines that pass the threshold!
#            #print(100*(contig_pass_coverage_threshold/(contig_pass_coverage_threshold+contig_fail_coverage_threshold)), "percent of", context, "cytosines pass coverage threshold!")
#            #print("Total methylated", context+"'s:", contig_cumulative_methylated)
#            #print("Total unmethylated", context+"'s:", contig_cumulative_unmethylated)
#            #print("Overall methylation percentage for passing", context, "cytosines is", 100*(contig_cumulative_methylated/(contig_cumulative_methylated+contig_cumulative_unmethylated)))
#        else:
#            #print("No data for this contig!")
    return contig_sample_methylation_stats_dict

def contig_to_feature(contig_sample_methylation_stats_dict, feature_annotation_dict, feature_type):
    contig_feature_cpg_coordinate_methylation_dict = {}
    contig_count = 0
    for i in contig_sample_methylation_stats_dict:
        contig_count +=1
        #print("Processing contig", contig_count, "Contig name:", i)
        coordinate_list = contig_sample_methylation_stats_dict[i][2]
        percent_list = contig_sample_methylation_stats_dict[i][3]
        methylated_list = contig_sample_methylation_stats_dict[i][4] # doesn't exist yet!
        unmethylated_list = contig_sample_methylation_stats_dict[i][5] # doesn't exist yet!
        feature_cpg_coordinate_dict = {}
        for j in range(len(feature_annotation_dict[i][0])): # start coordinate list
            start = int(feature_annotation_dict[i][0][j])
            end = int(feature_annotation_dict[i][1][j])
            if len(feature_annotation_dict[i]) == 3 and feature_type == "feature":
                strand = feature_annotation_dict[i][2][j]
            else:
                strand = "."
            coordinate = []
            percent = []
            methylated = []
            unmethylated = []
            #print("Start:", start, "End:", end, "Strand:", strand)
            for k in range(len(coordinate_list)):
                if coordinate_list[k] >= start and coordinate_list[k] <= end: # i.e. within feature
                    coordinate.append(coordinate_list[k])
                    methylated.append(methylated_list[k])
                    unmethylated.append(unmethylated_list[k])
                    percent.append(percent_list[k])
                    #print("coordinate:", coordinate_list[k], "percent:", percent_list[k], "methylated:", methylated_list[k], "unmethylated:", unmethylated_list[k])
            if strand == "+" or strand == "-":
                feature_cpg_coordinate_dict[i+"."+str(start)] = [start, end, strand, coordinate, percent, methylated, unmethylated]
            else: # strand == ".": # i.e. intergenic region
                if (end - (start-1))%2 == 0: #i.e. even number
                    fiveprime_feature_start = start
                    fiveprime_feature_end = int(start+((end-(start-1))/2))
                    threeprime_feature_start = int(start+(((end-(start-1))/2)+1))
                    threeprime_feature_end = end
                    fiveprime_feature_coordinate = []
                    threeprime_feature_coordinate = []
                    fiveprime_feature_methylated = []
                    threeprime_feature_methylated = []
                    fiveprime_feature_unmethylated = []
                    threeprime_feature_unmethylated = []
                    fiveprime_feature_percent = []
                    threeprime_feature_percent = []
                    for l in range(len(coordinate)):
                        if coordinate[l] >= fiveprime_feature_start and coordinate[l] <= fiveprime_feature_end:
                            fiveprime_feature_coordinate.append(coordinate[l])
                            fiveprime_feature_percent.append(percent[l])
                            fiveprime_feature_methylated.append(methylated[l])
                            fiveprime_feature_unmethylated.append(unmethylated[l])
                        elif coordinate[l] >= threeprime_feature_start and coordinate[l] <= threeprime_feature_end:
                            threeprime_feature_coordinate.append(coordinate[l])
                            threeprime_feature_percent.append(percent[l])
                            threeprime_feature_methylated.append(methylated[l])
                            threeprime_feature_unmethylated.append(unmethylated[l])
#                        else:
#                            print("Coordinate not found!")
                    feature_cpg_coordinate_dict[i+"."+str(start)] = \
[[fiveprime_feature_start, threeprime_feature_start], [fiveprime_feature_end, threeprime_feature_end], strand, [fiveprime_feature_coordinate, threeprime_feature_coordinate], \
[fiveprime_feature_percent, threeprime_feature_percent], [fiveprime_feature_methylated, threeprime_feature_methylated], [fiveprime_feature_unmethylated, threeprime_feature_unmethylated]]
                else: # odd number
                    fiveprime_feature_start = start
                    fiveprime_feature_end = int(start+(((end-(start-1))-1)/2))
                    threeprime_feature_start = int(start+((((end-(start-1))-1)/2)+2)) # plus two instead of 1!
                    threeprime_feature_end = end
                    fiveprime_feature_coordinate = []
                    threeprime_feature_coordinate = []
                    fiveprime_feature_methylated = []
                    threeprime_feature_methylated = []
                    fiveprime_feature_unmethylated = []
                    threeprime_feature_unmethylated = []
                    fiveprime_feature_percent = []
                    threeprime_feature_percent = []
                    for l in range(len(coordinate)):
                        if coordinate[l] >= fiveprime_feature_start and coordinate[l] <= fiveprime_feature_end:
                            fiveprime_feature_coordinate.append(coordinate[l])
                            fiveprime_feature_percent.append(percent[l])
                            fiveprime_feature_methylated.append(methylated[l])
                            fiveprime_feature_unmethylated.append(unmethylated[l])
                        elif coordinate[l] >= threeprime_feature_start and coordinate[l] <= threeprime_feature_end:
                            threeprime_feature_coordinate.append(coordinate[l])
                            threeprime_feature_percent.append(percent[l])
                            threeprime_feature_methylated.append(methylated[l])
                            threeprime_feature_unmethylated.append(unmethylated[l])
                        else:
#                            print("Coordinate not found!")
#                    feature_cpg_coordinate_dict[i+"."+str(start)] = \
[[fiveprime_feature_start, threeprime_feature_start], [fiveprime_feature_end, threeprime_feature_end], strand, [fiveprime_feature_coordinate, threeprime_feature_coordinate], \
[fiveprime_feature_percent, threeprime_feature_percent], [fiveprime_feature_methylated, threeprime_feature_methylated], [fiveprime_feature_unmethylated, threeprime_feature_unmethylated]]
        contig_feature_cpg_coordinate_methylation_dict[i] = feature_cpg_coordinate_dict
    return contig_feature_cpg_coordinate_methylation_dict

def add_upstream_downstream_regions(contig_gene_cpg_coordinate_methylation_dict, contig_intergene_cpg_coordinate_methylation_dict):
    contig_gene_and_flanking_sequences_cpg_coordinate_methylation_dict = {}
    for i in contig_gene_cpg_coordinate_methylation_dict:
        gene_and_flanking_sequences_cpg_coordinate_methylation_dict = {}
        for j in contig_gene_cpg_coordinate_methylation_dict[i]:
            gene_coordinate_list = contig_gene_cpg_coordinate_methylation_dict[i][j][3]
            gene_percent_list = contig_gene_cpg_coordinate_methylation_dict[i][j][4]
            gene_methylated_list = contig_gene_cpg_coordinate_methylation_dict[i][j][5]
            gene_unmethylated_list = contig_gene_cpg_coordinate_methylation_dict[i][j][6]
            combined_coordinate_list = [] # placeholder
            combined_percent_list = [] # placeholder
            combined_methylated_list = [] # placeholder
            combined_unmethylated_list = [] # placeholder
            for k in contig_intergene_cpg_coordinate_methylation_dict[i]:
                upstream_coordinate_list = [] # placeholder
                upstream_percent_list = [] # placeholder
                upstream_methylated_list = [] # placeholder
                upstream_unmethylated_list = [] # placeholder
                downstream_coordinate_list = [] # placeholder
                downstream_percent_list = [] # placeholder
                downstream_methylated_list = [] # placeholder
                downstream_unmethylated_list = [] # placeholder
                if (contig_intergene_cpg_coordinate_methylation_dict[i][k][1][1]) + 1 == contig_gene_cpg_coordinate_methylation_dict[i][j][0]: # if intergene three_prime_feature_end + 1 equals gene start (must be upstream)
                    upstream_coordinate_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][3][1]
                    upstream_percent_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][4][1]
                    upstream_methylated_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][5][1]
                    upstream_unmethylated_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][6][1]
                if (contig_intergene_cpg_coordinate_methylation_dict[i][k][0][0]) - 1 == contig_gene_cpg_coordinate_methylation_dict[i][j][1]: # if intergene five_prime_feature_start - 1 equals gene end (must be downstream)
                    downstream_coordinate_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][3][0]
                    downstream_percent_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][4][0]
                    downstream_methylated_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][5][0]
                    downstream_unmethylated_list = contig_intergene_cpg_coordinate_methylation_dict[i][k][6][0]
            combined_coordinate_list = upstream_coordinate_list + gene_coordinate_list + downstream_coordinate_list
            combined_percent_list = upstream_percent_list + gene_percent_list + downstream_percent_list
            combined_methylated_list = upstream_methylated_list + gene_methylated_list + downstream_methylated_list
            combined_unmethylated_list = upstream_unmethylated_list + gene_unmethylated_list + downstream_unmethylated_list
            gene_and_flanking_sequences_cpg_coordinate_methylation_dict[j] = \
[contig_gene_cpg_coordinate_methylation_dict[i][j][0], contig_gene_cpg_coordinate_methylation_dict[i][j][1], contig_gene_cpg_coordinate_methylation_dict[i][j][2], \
combined_coordinate_list, combined_percent_list, combined_methylated_list, combined_unmethylated_list]
        contig_gene_and_flanking_sequences_cpg_coordinate_methylation_dict[i] = gene_and_flanking_sequences_cpg_coordinate_methylation_dict
    return contig_gene_and_flanking_sequences_cpg_coordinate_methylation_dict

def gene_coordinate_to_position(contig_gene_cpg_coordinate_methylation_dict, zeropoint):
    contig_gene_cpg_position_methylation_dict = {}
    for i in contig_gene_cpg_coordinate_methylation_dict:
        gene_dict = {}
        for j in contig_gene_cpg_coordinate_methylation_dict[i]:
            #print(contig_gene_cpg_coordinate_methylation_dict[i][j][0], contig_gene_cpg_coordinate_methylation_dict[i][j][1], contig_gene_cpg_coordinate_methylation_dict[i][j][2])
            strand = contig_gene_cpg_coordinate_methylation_dict[i][j][2]
            start = contig_gene_cpg_coordinate_methylation_dict[i][j][0]
            end = contig_gene_cpg_coordinate_methylation_dict[i][j][1]
            if strand == "-":
                if len(contig_gene_cpg_coordinate_methylation_dict[i][j][3]) > 0:
                    first = contig_gene_cpg_coordinate_methylation_dict[i][j][3][0]
                    last = contig_gene_cpg_coordinate_methylation_dict[i][j][3][-1]
                    distance_of_start_from_first = start - first
                    distance_of_end_from_first = end - first
                    newstart = first + distance_of_start_from_first
                    newend = first + distance_of_end_from_first
                else:
                    continue
            if strand == "+":
                if zeropoint == "start":
                    zero = start
                elif zeropoint == "end":
                    zero = end
#                else:
#                    #print("Zeropoint error!")
            elif strand == "-":
                if zeropoint == "start":
                    zero = newstart
                elif zeropoint == "end":
                    zero = newend
#                else:
#                    #print("Zeropoint error!")
            cpg_position = []
            cpg_percent = []
            cpg_methylated = []
            cpg_unmethylated = []
            for k in range(len(contig_gene_cpg_coordinate_methylation_dict[i][j][3])): # scroll through coordinates for given gene
                if strand == "+":
                    position = contig_gene_cpg_coordinate_methylation_dict[i][j][3][k] - zero # this will mark the position 0
                    percent = contig_gene_cpg_coordinate_methylation_dict[i][j][4][k]
                    methylated = contig_gene_cpg_coordinate_methylation_dict[i][j][5][k]
                    unmethylated = contig_gene_cpg_coordinate_methylation_dict[i][j][6][k]
                    cpg_position.append(position)
                    cpg_percent.append(percent)
                    cpg_methylated.append(methylated)
                    cpg_unmethylated.append(unmethylated)
                elif strand == "-":
                    #print("position =", length, "- (", contig_gene_cpg_coordinate_methylation_dict[i][j][3][k], "-", start, ")")
                    distance_of_k_from_first = contig_gene_cpg_coordinate_methylation_dict[i][j][3][k] - first
                    flip = last - distance_of_k_from_first # flip coordinates around
                    position = flip - zero # now calculate position as above
                    percent = contig_gene_cpg_coordinate_methylation_dict[i][j][4][k]
                    methylated = contig_gene_cpg_coordinate_methylation_dict[i][j][5][k]
                    unmethylated = contig_gene_cpg_coordinate_methylation_dict[i][j][6][k]
                    cpg_position.append(position)
                    cpg_percent.append(percent)
                    cpg_methylated.append(methylated)
                    cpg_unmethylated.append(unmethylated)
#                else:
#                    #print("Error!")
            gene_dict[j] = [strand, start, end, cpg_position, cpg_percent, cpg_methylated, cpg_unmethylated]
        contig_gene_cpg_position_methylation_dict[i] = gene_dict
    return contig_gene_cpg_position_methylation_dict

def get_max_cpg_position(contig_gene_cpg_position_methylation_dict):
    contig_max_cpg_position = {}
    for i in contig_gene_cpg_position_methylation_dict:
        max_cpg_position = 0
        for j in contig_gene_cpg_position_methylation_dict[i]:
            if len(contig_gene_cpg_position_methylation_dict[i][j][3]) > 0:
                if max(contig_gene_cpg_position_methylation_dict[i][j][3]) > max_cpg_position:
                    max_cpg_position = max(contig_gene_cpg_position_methylation_dict[i][j][3])
        contig_max_cpg_position[i] = max_cpg_position
    return contig_max_cpg_position

def get_min_cpg_position(contig_gene_cpg_position_methylation_dict, contig_max_cpg_position):
    contig_min_cpg_position = {}
    for i in contig_gene_cpg_position_methylation_dict:
        min_cpg_position = contig_max_cpg_position[i]
        for j in contig_gene_cpg_position_methylation_dict[i]:
            if len(contig_gene_cpg_position_methylation_dict[i][j][3]) > 0:
                if min(contig_gene_cpg_position_methylation_dict[i][j][3]) < min_cpg_position:
                    min_cpg_position = min(contig_gene_cpg_position_methylation_dict[i][j][3])
        contig_min_cpg_position[i] = min_cpg_position
    return contig_min_cpg_position

def position_methylation_stat_tally(contig_min_position_dict, contig_max_position_dict, contig_feature_position_context_methylation_dict, stat):
    if stat == "percent":
        stat = 4
    if stat == "methylated_count":
        stat = 5
    if stat == "unmethylated_count":
        stat = 6
    statindex = stat
    position_stat_tally_dict = {}
    for i in range(min(contig_min_position_dict.values()), max(contig_max_position_dict.values())+1):
        position_stat_tally_dict[i] = []
    for i in contig_feature_position_context_methylation_dict:
        for j in contig_feature_position_context_methylation_dict[i]:
            if len(contig_feature_position_context_methylation_dict[i][j][3]) > 0:
                for k in range(len(contig_feature_position_context_methylation_dict[i][j][3])):
                    #print("position:", contig_feature_position_context_methylation_dict[i][j][3][k])
                    position = contig_feature_position_context_methylation_dict[i][j][3][k]
                    position_stat = contig_feature_position_context_methylation_dict[i][j][statindex][k]
                    position_stat_list = position_stat_tally_dict[position]
                    position_stat_list.append(position_stat)
                    position_stat_tally_dict[position] = position_stat_list
    return position_stat_tally_dict

def position_methylation_stats_to_list(percent_position_tally_dict, methylated_count_position_tally_dict, unmethylated_count_position_tally_dict):
    position_length = []
    position_position = []
    position_median_percent = []
    position_max_percent = []
    position_methylated_count = []
    position_unmethylated_count = []
    for i in percent_position_tally_dict:
        position_length.append(len(percent_position_tally_dict[i]))
        position_position.append(i)
        if len(percent_position_tally_dict[i]) > 0:
            position_median_percent.append(statistics.median(percent_position_tally_dict[i]))
            if len(percent_position_tally_dict) > 1:
                position_max_percent.append(max(percent_position_tally_dict[i]))
            else:
                position_max_percent.append(percent_position_tally_dict[i][0])
        else:
            position_median_percent.append(-1)
            position_max_percent.append(-1)
    for i in methylated_count_position_tally_dict:
        if len(methylated_count_position_tally_dict[i]) > 0:
            methylated_count = 0
            for j in methylated_count_position_tally_dict[i]:
                methylated_count += j
            position_methylated_count.append(methylated_count)
        else:
            position_methylated_count.append(-1)
    for i in unmethylated_count_position_tally_dict:
        if len(unmethylated_count_position_tally_dict[i]) > 0:
            unmethylated_count = 0
            for j in unmethylated_count_position_tally_dict[i]:
                unmethylated_count += j
            position_unmethylated_count.append(unmethylated_count)
        else:
            position_unmethylated_count.append(-1)
    stats_list = [position_position, position_length, position_median_percent, position_max_percent, position_methylated_count, position_unmethylated_count]
    return stats_list

def get_feature_methylation_counts(contig_feature_position_methylation_dict, start_position, end_position):
    # record methylation data between selected positions for each gene, i.e. average methylation
    contig_feature_selected_positions_stats_dict = {}
    for i in contig_feature_position_methylation_dict:
        feature_selected_positions_stats_dict = {}
        feature_selected_methylated_count = 0
        feature_selected_unmethylated_count = 0
        for j in contig_feature_position_methylation_dict[i]:
            feature_selected_positions_methylated = []
            feature_selected_positions_unmethylated = []
            for k in range(len(contig_feature_position_methylation_dict[i][j][3])):
                if contig_feature_position_methylation_dict[i][j][3][k] <= end_position and contig_feature_position_methylation_dict[i][j][3][k] >= start_position:
                    feature_selected_positions_methylated.append(contig_feature_position_methylation_dict[i][j][5][k])
                    feature_selected_positions_unmethylated.append(contig_feature_position_methylation_dict[i][j][6][k])
            for m in feature_selected_positions_methylated:
                feature_selected_methylated_count += m
            for u in feature_selected_positions_unmethylated:
                feature_selected_unmethylated_count += u
            feature_selected_positions_stats_dict[j] = [feature_selected_methylated_count, feature_selected_unmethylated_count]
        contig_feature_selected_positions_stats_dict[i] = feature_selected_positions_stats_dict
    return contig_feature_selected_positions_stats_dict

def feature_methylation_counts_to_percent_list(contig_feature_selected_positions_stats_dict):
    selected_feature_position_methylation_list = []
    for i in contig_feature_selected_positions_stats_dict:
        for j in contig_feature_selected_positions_stats_dict[i]:
            if contig_feature_selected_positions_stats_dict[i][j][1] > 0:
                selected_feature_position_methylation_list.append(100*(contig_feature_selected_positions_stats_dict[i][j][0]/(contig_feature_selected_positions_stats_dict[i][j][0]+contig_feature_selected_positions_stats_dict[i][j][1])))
    return selected_feature_position_methylation_list
