####### IMPORT MODULES ########

import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mt_functions import *

### CREATE A CONTIG SEQUENCE DICTIONARY ###

#contig.genome_version
#contig.name = # name
#contig.size = # numerical value
#contig.gc # gc content

path_to_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest"
cx_report = "adultCombined.CX_report.txt"
sample_CX_data_file = path_to_data+"/"+cx_report
fasta = "GCF_002156985.1_Harm_1.0_genomic.fa"
gff = "GCF_002156985.1_Harm_1.0_genomic.gff"

contig_list = ["NC_014668.1","NW_018396379.1","NW_018396386.1","NW_018396377.1","NW_018396384.1"]

contig_sequence_dict = import_contig_sequence(contig_list)

print(len(contig_sequence_dict), "contigs to be analysed!")

for i in contig_sequence_dict:
    print(i,"\t", contig_sequence_dict[i][0:101]) # print first 100 nucleotides of each contig

gene_annotation = import_gff_data(contig_list, "gene")

intergene_annotation = import_intergene_coordinates(gene_annotation) # find non genic coordinates in gff format (1-based)

adult_combined_methylation_dict = import_methylation_data(contig_list, sample_CX_data_file)

cpg_coordinates = get_coordinates(contig_sequence_dict, "CG")
# check that it works!
# contig_sequence_dict["NC_014668.1"][cpg_coordinates["NC_014668.1"][0][2]:cpg_coordinates["NC_014668.1"][1][2]]

# extract gene stats (i.e. size, cpg count, cpg percent, cpg coordinates)
gene_cpg_count_dict = {}
gene_cpg_coords_dict = {}
gene_size_dict = {}
gene_cpg_percent_dict = {}

for i in contig_sequence_dict:
    number_of_cpgs = []
    cpg_coords = []
    size_of_gene = []
    gene_percent_cpg = []
    for j in range(len(gene_annotation[i][0])):
        gene_start_coordinate = int(gene_annotation[i][0][j])
        gene_end_coordinate = int(gene_annotation[i][1][j])
        gene_cpg_count = 0
        gene_cpg_coords = []
        for k in range(len(cpg_coordinates[i][0])):
            if int(cpg_coordinates[i][0][k]) <= gene_end_coordinate:
                if int(cpg_coordinates[i][0][k]) >= gene_start_coordinate:
                    gene_cpg_count += 1
                    gene_cpg_coords.append(cpg_coordinates[i][0][k])
                    print(i,"\t",cpg_coordinates[i][0][k],"\t",contig_sequence_dict[i][cpg_coordinates[i][0][k]:cpg_coordinates[i][1][k]])
        number_of_cpgs.append(gene_cpg_count)
        cpg_coords.append(gene_cpg_coords)
        size_of_gene.append(gene_end_coordinate-gene_start_coordinate)
        gene_percent_cpg.append(100*(gene_cpg_count/(gene_end_coordinate-gene_start_coordinate)))
        print("gene_cpg_count:", gene_cpg_count)
    gene_cpg_count_dict[i] = number_of_cpgs
    gene_cpg_coords_dict[i] = cpg_coords
    gene_size_dict[i] = size_of_gene
    gene_cpg_percent_dict[i] = gene_percent_cpg

#contig.annotation = # gene annotation, CpG islands etc

### CREATE A FEATURE OBJECT (INHERIT FROM CONTIG) ###
#feature.type
#feature.start
#feature.end

### CREATE A GENE OBJECT (INHERIT FROM FEATURE) ###
#gene.start = # start
#gene.end = # end
#gene.tss =  # transcription start site
#gene.strand = # plus or minus
#gene.cds = # coding sequence
#gene.annotation = # exons, introns

### CREATE A SAMPLE OBJECT ###
#sample.methylation_coordinates
#sample.sequence_map
#sample.coverage

#### CALCULATE METHYLATION STATS FOR EACH CONTIG ########

contig_gene_methylation_dict = {}

for i in contig_list:
    cumulative_methylated = 0
    cumulative_unmethylated = 0
    for j in gene_cpg_coords_dict[i]:
        for k in j:
            coordinate = adult_combined_methylation_dict[i][0][adult_combined_methylation_dict[i][0].index(k)]
            methylated = adult_combined_methylation_dict[i][2][adult_combined_methylation_dict[i][0].index(k)]
            cumulative_methylated += adult_combined_methylation_dict[i][2][adult_combined_methylation_dict[i][0].index(k)]
            unmethylated = adult_combined_methylation_dict[i][3][adult_combined_methylation_dict[i][0].index(k)]
            cumulative_unmethylated += adult_combined_methylation_dict[i][3][adult_combined_methylation_dict[i][0].index(k)]
            total = methylated+unmethylated
            if total > 0:
                percent = 100*(methylated/total)
            else:
                percent = "No data!"
            context = adult_combined_methylation_dict[i][4][adult_combined_methylation_dict[i][0].index(k)]
            print("coordinate:", coordinate, "methylated:", methylated, "total:", total, "percent:", percent, "context:", context)
    totals_data = [cumulative_methylated, cumulative_unmethylated]
    #print("cumulative percent calculation 1:", 100*(cumulative_methylated/(cumulative_methylated+cumulative_unmethylated)))
    contig_gene_methylation_dict[i] = totals_data

#########################################################

####### set upstream, downstream, 5prime and 3prime lengths #######

#upstream_length = 1000
#downstream_length = 1000
#fiveprime_length = 1000
#threeprime_length = 1000

####### set total methylation tallys #######

#total_methylated_CpG_tally = {}
#total_unmethylated_CpG_tally = {}
#total_upstream_methylated_CpG_tally = {}
#total_upstream_unmethylated_CpG_tally = {}

#for i in range(1,threeprime_length+1):
#    total_unmethylated_CpG_tally[i] = 0
#    total_methylated_CpG_tally[i] = 0

#for i in range(1,upstream_length+1):
#    total_upstream_methylated_CpG_tally[i] = 0
#    total_upstream_unmethylated_CpG_tally[i] = 0

# initialize gene number tally to zero

#gene_number = 0

####### loop through contigs #######

#for entry in contig_list: # add indices to restrict to particular contigs in list
#
#    print("Contig name:", entry)
#    contig=open("test_contig_sequence/"+entry+"_sequence.txt", "r") # open the file containing contig sequence
#    contig_sequence=contig.read() # and save to variable 'contig_sequence'
#    gene_start=open("test_contig_beds_annotations/"+entry+"_gene_start.txt", 'r') # open file that contains gene start information for specific contig
#    gene_end=open("test_contig_beds_annotations/"+entry+"_gene_end.txt", 'r') # and the file that contains gene end information
#    gene_strand=open("test_contig_beds_annotations/"+entry+"_gene_strand.txt", 'r') # and the file that contains strand information

#    gene_start_coordinates=[] # create an empty list to store the gene_start coordinates
#    gene_end_coordinates=[] # and to store the gene_end_coordinates
#    gene_strand_designation=[] # and to store the gene strand

#    gene_5prime_CpG_positions={} # make a dictionary to store the CpG coordinates occurring at the start (~ 500 bp) of each gene
#    gene_3prime_CpG_positions={} # and one to store CpG coordinates occurring at the end (~500 bp) of each gene
#    upstream_CpG_positions = {} ### initialize upstream CpG dictionary

####### CONTIG DATA EXTRACTION #######

#    contig_CpG_position, contig_GpC_position, contig_CHG_position, contig_CHH_position = find_methylation_sites(contig_sequence)

#    print("CpG count for contig:", len(contig_CpG_position))
#    print("GC content for contig:",((len(contig_CpG_position)*2)+(len(contig_GpC_position)*2))/len(contig_sequence))

#    for i in gene_start:
#        i.strip()
#        gene_start_coordinates.append(int(i)) # store gene_start coordinates
#    for i in gene_end:
#        i.strip()
#        gene_end_coordinates.append(int(i)) # store gene_end coordinates
#    for i in gene_strand:
#        designation=i.strip()
#        gene_strand_designation.append(designation) # store gene_strand

#    print("Number of genes in contig:", len(gene_start_coordinates), "\t", gene_start_coordinates[0:6], "\t", gene_end_coordinates[0:6], "\t", gene_strand_designation[0:6])

#    contig.close()
#    gene_start.close()
#    gene_end.close()
#    gene_strand.close()

####### GENE DATA EXTRACTION #######

#    for i in range(len(gene_start_coordinates)): # loop through each gene in contig
#        gene_number += 1
#        print("Processing gene", gene_number)
        #print("gene start coordinate:", gene_start_coordinates[i], "\t", "gene end coordinate:", gene_end_coordinates[i], "\t", "strand:", gene_strand_designation[i]) # print the start and end coordiantes for this gene
        #print(contig_sequence[gene_start_coordinates[i]:gene_end_coordinates[i]]) # print the gene sequence
#        if gene_strand_designation[i]=="+": # if the gene is on the plus strand
#            gene_sequence = contig_sequence[gene_start_coordinates[i]:gene_end_coordinates[i]] # save the sequence in the orientation found in the contig_sequence variable
#            if gene_start_coordinates[i]-upstream_length >= 0:
#                upstream_sequence = contig_sequence[gene_start_coordinates[i]-upstream_length:gene_start_coordinates[i]]
#            else:
#                upstream_sequence = contig_sequence[0:gene_start_coordinates[i]]
#        else:
#            gene_sequence = Seq(contig_sequence[gene_start_coordinates[i]:gene_end_coordinates[i]], generic_dna).reverse_complement() # or save the reverse complement sequence if the gene is found on the minus strand
#            if gene_end_coordinates[i]+upstream_length < len(contig_sequence)+1:
#                upstream_sequence = Seq(contig_sequence[gene_end_coordinates[i]:gene_end_coordinates[i]+upstream_length], generic_dna).reverse_complement()
#            else:
#                upstream_sequence = Seq(contig_sequence[gene_end_coordinates[i]:len(contig_sequence)+1], generic_dna).reverse_complement()


        ####### Find the CpG, GpC positions in the gene sequence #######

#        CpG_position, GpC_position, CHG_position, CHH_position = find_methylation_sites(gene_sequence)

        ####### Find the CpG, GpC positions in the upstream sequence #######

#        upstream_CpG_position, upstream_GpC_position, upstream_CHG_position, upstream_CHH_position = find_methylation_sites(upstream_sequence)

        ####### Find the CpG, GpC positions in the 5prime and 3prime sequences #######

#        if len(gene_sequence) > fiveprime_length+threeprime_length:
#            threeprime_CpG_position = find_methylation_sites(gene_sequence[0:threeprime_length+1])
#        else:
#            threeprime_CpG_position = []
#        if len(gene_sequence) > fiveprime_length+threeprime_length:
#            fiveprime_CpG_position = find_methylation_sites(gene_sequence[len(gene_sequence)-fiveprime_length:len(gene_sequence)+1])
#        else:
#            fiveprime_CpG_position = []

        ####### Add the position lists to relevant dictionary #######

#        gene_5prime_CpG_positions[gene_start_coordinates[i]]=fiveprime_CpG_position # add start CpG list to a dictionary with gene start coordinate as key
#        gene_3prime_CpG_positions[gene_end_coordinates[i]]=threeprime_CpG_position # same as above for end CpG list

#        upstream_CpG_positions[gene_start_coordinates[i]]=upstream_CpG_position

####### Print line in summary table #######

####### DATASET SELECTION #######

### Import the methylation data for the contig ###

#    zero_cov_coordinates = open("input_files/"+entry+"_coordinate.txt","r")
#    zero_cov_methylated = open("input_files/"+entry+"_methylated.txt","r")
#    zero_cov_unmethylated = open("input_files/"+entry+"_unmethylated.txt","r")

#    coordinate_list=[]
#    methylated_list=[]
#    unmethylated_list=[]

#    for i in zero_cov_coordinates:
#        coordinate = int(i.strip())
#        coordinate_list.append(coordinate)

#    for i in zero_cov_methylated:
#        methylated = int(i.strip())
#        methylated_list.append(methylated)

#    for i in zero_cov_unmethylated:
#        unmethylated = int(i.strip())
#        unmethylated_list.append(unmethylated)

#    zero_cov_coordinates.close()
#    zero_cov_methylated.close()
#    zero_cov_unmethylated.close()

#    for i in range(0,10):
#        print(coordinate_list[i], contig_sequence[i:i+2])

### Split contig methylation data into genes (i.e. gene start coordinate as dictionary key) ###

#    zero_cov_data={}

#    for i in range(len(coordinate_list)):
#        zero_cov_data[coordinate_list[i]] = [methylated_list[i], unmethylated_list[i]]

#    for i in gene_start_coordinates: # loop through the gene start coordinates for each gene in the contig

#        methylated_CpG_tally = {}
#        unmethylated_CpG_tally = {}
#        no_CpG_data=[]

#        upstream_methylated_CpG_tally = {}
#        upstream_unmethylated_CpG_tally = {}
#        upstream_no_CpG_data = []

        # initialize methylated and unmethylated tally dictionaries to start at 0 for each position

#        for j in range(1,threeprime_length+1):
#            unmethylated_CpG_tally[j] = 0
#        for j in range(1,threeprime_length+1):
#            methylated_CpG_tally[j] = 0
#        for j in range(1,upstream_length+1):
#            upstream_unmethylated_CpG_tally[j] = 0
#            upstream_methylated_CpG_tally[j] = 0

#        gene_start_coordinate = i # assign gene start coordinate value to variable 'gene_start_coordinate'
#        coordinate_CpG_positions = gene_5prime_CpG_positions[i]

#        upstream_coordinate_CpG_positions = upstream_CpG_positions[i]

#        for j in coordinate_CpG_positions: # loop through the list of CpG positions
#            CpG_coordinate = j+gene_start_coordinate # convert position into coordinate
#            print("coordinate", CpG_coordinate, ":", gene_sequence[j:j+2]) ## check we're hitting CpG's!!!
#            if CpG_coordinate in coordinate_list:
#                methylated_CpG_tally[j] += zero_cov_data[CpG_coordinate][0] # add coordinate methylation value to the relevent position tally
#                unmethylated_CpG_tally[j] += zero_cov_data[CpG_coordinate][1] # add coordinate unmethylation value to the relevant position tally
#            else:
#                no_CpG_data.append(j)
#        for j in upstream_coordinate_CpG_positions:
#            upstream_CpG_coordinate = (gene_start_coordinate-upstream_length)+j #### CHECK THIS!!! ####
#            if CpG_coordinate in coordinate_list:
#                upstream_methylated_CpG_tally[j] += zero_cov_data[CpG_coordinate][0]
#                upstream_unmethylated_CpG_tally[j] += zero_cov_data[CpG_coordinate][1]
#            else:
#                upstream_no_CpG_data.append(j)

        ### Retrieve per gene methylation data metrics ###

#        percent_CpG_methylation_string = ""
#        CpG_coverage_string = ""

#        for j in range(1,threeprime_length+1):
#            position_coverage = methylated_CpG_tally[j]+unmethylated_CpG_tally[j]
#            if methylated_CpG_tally[j] > 0:
#                percent_methylation = 100*(methylated_CpG_tally[j]/position_coverage)
#            else:
#                percent_methylation = 0
#            percent_CpG_methylation_string = percent_CpG_methylation_string+str(percent_methylation)+"\t"
#            CpG_coverage_string = CpG_coverage_string+str(position_coverage)+"\t"

#        methylation_output.write(entry+"_"+str(i)+"\t"+percent_CpG_methylation_string+"\n")
#        coverage_output.write(entry+"_"+str(i)+"\n"+CpG_coverage_string+"\n")

        ### Retrieve per gene dictionary data ###

#        for i in methylated_CpG_tally:
#            total_methylated_CpG_tally[i] += methylated_CpG_tally[i]

#        for i in unmethylated_CpG_tally:
#            total_unmethylated_CpG_tally[i] += unmethylated_CpG_tally[i]

        ### Retrieve upstream per gene dictionary data ###

#        for i in upstream_unmethylated_CpG_tally:
#            total_upstream_unmethylated_CpG_tally[i] += upstream_unmethylated_CpG_tally[i]

#        for i in upstream_methylated_CpG_tally:
#            total_upstream_methylated_CpG_tally[i] += upstream_methylated_CpG_tally[i]

#print("total methylated CpG tally:", total_methylated_CpG_tally)
#print("total unmethylated CpG tally:", total_unmethylated_CpG_tally)
#print("total upstream methylated CpG tally:", total_upstream_methylated_CpG_tally)
#print("total upstream unmethylated CpG tally:", total_upstream_unmethylated_CpG_tally)
