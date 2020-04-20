####### IMPORT MODULES ########

import matplotlib as mpl
import numpy as np
mpl.use('agg')
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mt_functions import *
from plot_functions import *
import statistics

### DEFINE PATHS TO DATA FILES ###

path_to_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest"
fasta = "GCF_002156985.1_Harm_1.0_genomic.fa"
gff = "GCF_002156985.1_Harm_1.0_genomic.gff"
P4_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-P4/Jack-methylation-P4_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
P5_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-P5/Jack-methylation-P5_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
P6_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-P6/Jack-methylation-P6_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
pupa_combined_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/methylationExtractPupaCombine/pupaCombined.CX_report.txt"
adult_combined_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/methylationExtractAdultCombine/adultCombined.CX_report.txt"
P4_SSSI_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-SSS1-P4/Jack-methylation-repli-g-SSS1-P4_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
A7_SSSI_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-A7/Jack-methylation-repli-g-A7_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
P4_replig_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-P4/Jack-methylation-repli-g-P4_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
A7_replig_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-A7/Jack-methylation-repli-g-A7_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"

### CREATE A CONTIG SEQUENCE DICTIONARY ###

#contig_list = ["NC_014668.1","NW_018396379.1","NW_018396386.1","NW_018396377.1","NW_018396384.1"]
contig_list = import_contig_list("/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest/contigs.txt")
contig_sequence_dict = import_contig_sequence(contig_list)

print(len(contig_sequence_dict), "contigs to be analysed!")

for i in contig_sequence_dict:
    print(i,"\t", contig_sequence_dict[i][0:101]) # print first 100 nucleotides of each contig

### GET GENE ANNOTATION DATA ###

gene_annotation = import_gff_data(contig_list, "gene")

### GET NON-GENE COORDINATES ###

intergene_annotation = import_intergene_coordinates(gene_annotation, contig_sequence_dict) # find non genic coordinates in gff format (1-based)

### IMPORT METHYLATION DATA ###

P4_methylation_dict = import_methylation_data(contig_list, P4_CX_data)
P5_methylation_dict = import_methylation_data(contig_list, P5_CX_data)
P6_methylation_dict = import_methylation_data(contig_list, P6_CX_data)
pupa_combined_methylation_dict = import_methylation_data(contig_list, pupa_combined_CX_data)
adult_combined_methylation_dict = import_methylation_data(contig_list, adult_combined_CX_data)
P4_SSSI_methylation_dict = import_methylation_data(contig_list, P4_SSSI_CX_data)
A7_SSSI_methylation_dict = import_methylation_data(contig_list, A7_SSSI_CX_data)
P4_replig_methylation_dict = import_methylation_dict(contig_list, P4_replig_CX_data)
A7_replig_methylation_dict = import_methylation_data(contig_list, A7_replig_CX_data)

### MARK CPG COORDINATES IN IMPORTED CONTIGS ###

cpg_coordinates = get_coordinates(contig_sequence_dict, "CG")
# check that it works!
# contig_sequence_dict["NC_014668.1"][cpg_coordinates["NC_014668.1"][0][2]:cpg_coordinates["NC_014668.1"][1][2]]

### EXTRACT STATS FOR GENIC AND INTERGENIC REGIONS ###

# extract stats (i.e. feature size, cpg count, cpg percent, cpg coordinates)
contig_gene_stats_dict = get_feature_stats(cpg_coordinates, gene_annotation, contig_sequence_dict)
contig_intergene_stats_dict = get_feature_stats(cpg_coordinates, intergene_annotation, contig_sequence_dict)

### ANALYSE SAMPLES ###

pupa_contig_cg_methylation_stats_dict = get_sample_methylation_stats(pupa_combined_methylation_dict, 20, "CG")
pupa_contig_chg_methylation_stats_dict = get_sample_methylation_stats(pupa_combined_methylation_dict, 20, "CHG")
pupa_contig_chh_methylation_stats_dict = get_sample_methylation_stats(pupa_combined_methylation_dict, 20, "CHH")
adult_contig_cg_methylation_stats_dict = get_sample_methylation_stats(adult_combined_methylation_dict, 20, "CG")
adult_contig_chg_methylation_stats_dict = get_sample_methylation_stats(adult_combined_methylation_dict, 20, "CHG")
adult_contig_chh_methylation_stats_dict = get_sample_methylation_stats(adult_combined_methylation_dict, 20, "CHH")

# group data by gene and intergene (i.e. feature) coordinates
pupa_gene_cg_methylation = contig_to_feature(pupa_contig_cg_methylation_stats_dict, gene_annotation)
pupa_intergene_cg_methylation = contig_to_feature(pupa_contig_cg_methylation_stats_dict, intergene_annotation)

pupa_gene_chg_methylation = contig_to_feature(pupa_contig_chg_methylation_stats_dict, gene_annotation)
pupa_intergene_chg_methylation = contig_to_feature(pupa_contig_chg_methylation_stats_dict, intergene_annotation)

pupa_gene_chh_methylation = contig_to_feature(pupa_contig_chh_methylation_stats_dict, gene_annotation)
pupa_intergene_chh_methylation = contig_to_feature(pupa_contig_chh_methylation_stats_dict, intergene_annotation)

adult_gene_cg_methylation = contig_to_feature(adult_contig_cg_methylation_stats_dict, gene_annotation)
adult_intergene_cg_methylation = contig_to_feature(adult_contig_cg_methylation_stats_dict, intergene_annotation)

adult_gene_chg_methylation = contig_to_feature(adult_contig_chg_methylation_stats_dict, gene_annotation)
adult_intergene_chg_methylation = contig_to_feature(adult_contig_chg_methylation_stats_dict, intergene_annotation)

adult_gene_chh_methylation = contig_to_feature(adult_contig_chh_methylation_stats_dict, gene_annotation)
adult_intergene_chh_methylation = contig_to_feature(adult_contig_chh_methylation_stats_dict, intergene_annotation)

# add flanking regions to gene coordinate data (upstream and downstream)
pupa_gene_plus_flanking_cg_methylation = add_upstream_downstream_regions(pupa_gene_cg_methylation, pupa_intergene_cg_methylation)
pupa_gene_plus_flanking_chg_methylation = add_upstream_downstream_regions(pupa_gene_chg_methylation, pupa_intergene_chg_methylation)
pupa_gene_plus_flanking_chh_methylation = add_upstream_downstream_regions(pupa_gene_chh_methylation, pupa_intergene_chh_methylation)

adult_gene_plus_flanking_cg_methylation = add_upstream_downstream_regions(adult_gene_cg_methylation, adult_intergene_cg_methylation)
adult_gene_plus_flanking_chg_methylation = add_upstream_downstream_regions(adult_gene_chg_methylation, adult_intergene_chg_methylation)
adult_gene_plus_flanking_chh_methylation = add_upstream_downstream_regions(adult_gene_chh_methylation, adult_intergene_chh_methylation)

# check upstream, gene and downstream cpg count stats
#upstream = []
#downstream = []
#gene = []
#for i in adult_gene_plus_flanking_methylation:
#    for j in adult_gene_plus_flanking_methylation[i]:
#        for k in range(len(adult_gene_plus_flanking_methylation[i][j][3])):
#            if adult_gene_plus_flanking_methylation[i][j][3][k] < adult_gene_plus_flanking_methylation[i][j][0]:
#                upstream.append(adult_gene_plus_flanking_methylation[i][j][3][k])
#            elif adult_gene_plus_flanking_methylation[i][j][3][k] > adult_gene_plus_flanking_methylation[i][j][1]:
#                downstream.append(adult_gene_plus_flanking_methylation[i][j][3][k])
#            else:
#                gene.append(adult_gene_plus_flanking_methylation[i][j][3][k])
#
#len(upstream)
#len(gene)
#len(downstream)

# convert coordinates to positions
pupa_gene_plus_flanking_cg_methylation_position = gene_coordinate_to_position(pupa_gene_plus_flanking_cg_methylation, "start")
pupa_gene_plus_flanking_cg_methylation_position_end = gene_coordinate_to_position(pupa_gene_plus_flanking_cg_methylation, "end")
pupa_gene_plus_flanking_chg_methylation_position = gene_coordinate_to_position(pupa_gene_plus_flanking_chg_methylation, "start")
pupa_gene_plus_flanking_chg_methylation_position_end = gene_coordinate_to_position(pupa_gene_plus_flanking_chg_methylation, "end")
pupa_gene_plus_flanking_chh_methylation_position = gene_coordinate_to_position(pupa_gene_plus_flanking_chh_methylation, "start")
pupa_gene_plus_flanking_chh_methylation_position_end = gene_coordinate_to_position(pupa_gene_plus_flanking_chh_methylation, "end")

adult_gene_plus_flanking_cg_methylation_position = gene_coordinate_to_position(adult_gene_plus_flanking_cg_methylation, "start")
adult_gene_plus_flanking_cg_methylation_position_end = gene_coordinate_to_position(adult_gene_plus_flanking_cg_methylation, "end")
adult_gene_plus_flanking_chg_methylation_position = gene_coordinate_to_position(adult_gene_plus_flanking_chg_methylation, "start")
adult_gene_plus_flanking_chg_methylation_position_end = gene_coordinate_to_position(adult_gene_plus_flanking_chg_methylation, "end")
adult_gene_plus_flanking_chh_methylation_position = gene_coordinate_to_position(adult_gene_plus_flanking_chh_methylation, "start")
adult_gene_plus_flanking_chh_methylation_position_end = gene_coordinate_to_position(adult_gene_plus_flanking_chh_methylation, "end")

# find max and min positions
pupa_max_cpg_position = get_max_cpg_position(pupa_gene_plus_flanking_cg_methylation_position)
pupa_max_cpg_position_end = get_max_cpg_position(pupa_gene_plus_flanking_cg_methylation_position_end)
pupa_max_chg_position = get_max_cpg_position(pupa_gene_plus_flanking_chg_methylation_position)
pupa_max_chg_position_end = get_max_cpg_position(pupa_gene_plus_flanking_chg_methylation_position_end)
pupa_max_chh_position = get_max_cpg_position(pupa_gene_plus_flanking_chh_methylation_position)
pupa_max_chh_position_end = get_max_cpg_position(pupa_gene_plus_flanking_chh_methylation_position_end)

adult_max_cpg_position = get_max_cpg_position(adult_gene_plus_flanking_cg_methylation_position)
adult_max_cpg_position_end = get_max_cpg_position(adult_gene_plus_flanking_cg_methylation_position_end)
adult_max_chg_position = get_max_cpg_position(adult_gene_plus_flanking_chg_methylation_position)
adult_max_chg_position_end = get_max_cpg_position(adult_gene_plus_flanking_chg_methylation_position_end)
adult_max_chh_position = get_max_cpg_position(adult_gene_plus_flanking_chh_methylation_position)
adult_max_chh_position_end = get_max_cpg_position(adult_gene_plus_flanking_chh_methylation_position_end)

pupa_min_cpg_position = get_min_cpg_position(pupa_gene_plus_flanking_cg_methylation_position, pupa_max_cpg_position)
pupa_min_cpg_position_end = get_min_cpg_position(pupa_gene_plus_flanking_cg_methylation_position_end, pupa_max_cpg_position_end)
pupa_min_chg_position = get_min_cpg_position(pupa_gene_plus_flanking_chg_methylation_position, pupa_max_chg_position)
pupa_min_chg_position_end = get_min_cpg_position(pupa_gene_plus_flanking_chg_methylation_position_end, pupa_max_chg_position_end)
pupa_min_chh_position = get_min_cpg_position(pupa_gene_plus_flanking_chh_methylation_position, pupa_max_chh_position)
pupa_min_chh_position_end = get_min_cpg_position(pupa_gene_plus_flanking_chh_methylation_position_end, pupa_max_chh_position_end)

adult_min_cpg_position = get_min_cpg_position(adult_gene_plus_flanking_cg_methylation_position, adult_max_cpg_position)
adult_min_cpg_position_end = get_min_cpg_position(adult_gene_plus_flanking_cg_methylation_position_end, adult_max_cpg_position_end)
adult_min_chg_position = get_min_cpg_position(adult_gene_plus_flanking_chg_methylation_position, adult_max_chg_position)
adult_min_chg_position_end = get_min_cpg_position(adult_gene_plus_flanking_chg_methylation_position_end, adult_max_chg_position_end)
adult_min_chh_position = get_min_cpg_position(adult_gene_plus_flanking_chh_methylation_position, adult_max_chh_position)
adult_min_chh_position_end = get_min_cpg_position(adult_gene_plus_flanking_chh_methylation_position_end, adult_max_chh_position_end)

# create percent cpg dictionary entry for each position between min and max
pupa_percent_cg_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_cpg_position, pupa_max_cpg_position, pupa_gene_plus_flanking_cg_methylation_position, "percent")
pupa_methylated_count_cg_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_cpg_position, pupa_max_cpg_position, pupa_gene_plus_flanking_cg_methylation_position, "methylated_count")
pupa_unmethylated_count_cg_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_cpg_position, pupa_max_cpg_position, pupa_gene_plus_flanking_cg_methylation_position, "unmethylated_count")
pupa_percent_cg_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_cpg_position_end, pupa_max_cpg_position_end, pupa_gene_plus_flanking_cg_methylation_position_end, "percent")
pupa_methylated_count_cg_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_cpg_position_end, pupa_max_cpg_position_end, pupa_gene_plus_flanking_cg_methylation_position_end, "methylated_count")
pupa_unmethylated_count_cg_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_cpg_position_end, pupa_max_cpg_position_end, pupa_gene_plus_flanking_cg_methylation_position_end, "unmethylated_count")

pupa_percent_chg_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_chg_position, pupa_max_chg_position, pupa_gene_plus_flanking_chg_methylation_position, "percent")
pupa_methylated_count_chg_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_chg_position, pupa_max_chg_position, pupa_gene_plus_flanking_chg_methylation_position, "methylated_count")
pupa_unmethylated_count_chg_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_chg_position, pupa_max_chg_position, pupa_gene_plus_flanking_chg_methylation_position, "unmethylated_count")
pupa_percent_chg_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_chg_position_end, pupa_max_chg_position_end, pupa_gene_plus_flanking_chg_methylation_position_end, "percent")
pupa_methylated_count_chg_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_chg_position_end, pupa_max_chg_position_end, pupa_gene_plus_flanking_chg_methylation_position_end, "methylated_count")
pupa_unmethylated_count_chg_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_chg_position_end, pupa_max_chg_position_end, pupa_gene_plus_flanking_chg_methylation_position_end, "unmethylated_count")

pupa_percent_chh_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_chh_position, pupa_max_chh_position, pupa_gene_plus_flanking_chh_methylation_position, "percent")
pupa_methylated_count_chh_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_chh_position, pupa_max_chh_position, pupa_gene_plus_flanking_chh_methylation_position, "methylated_count")
pupa_unmethylated_count_chh_methylation_position_tally_dict = position_methylation_stat_tally(pupa_min_chh_position, pupa_max_chh_position, pupa_gene_plus_flanking_chh_methylation_position, "unmethylated_count")
pupa_percent_chh_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_chh_position_end, pupa_max_chh_position_end, pupa_gene_plus_flanking_chh_methylation_position_end, "percent")
pupa_methylated_count_chh_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_chh_position_end, pupa_max_chh_position_end, pupa_gene_plus_flanking_chh_methylation_position_end, "methylated_count")
pupa_unmethylated_count_chh_methylation_position_end_tally_dict = position_methylation_stat_tally(pupa_min_chh_position_end, pupa_max_chh_position_end, pupa_gene_plus_flanking_chh_methylation_position_end, "unmethylated_count")

adult_percent_cg_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_cpg_position, adult_max_cpg_position, adult_gene_plus_flanking_cg_methylation_position, "percent")
adult_methylated_cg_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_cpg_position, adult_max_cpg_position, adult_gene_plus_flanking_cg_methylation_position, "methylated_count")
adult_unmethylated_cg_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_cpg_position, adult_max_cpg_position, adult_gene_plus_flanking_cg_methylation_position, "unmethylated_count")
adult_percent_cg_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_cpg_position_end, adult_max_cpg_position_end, adult_gene_plus_flanking_cg_methylation_position_end, "percent")
adult_methylated_count_cg_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_cpg_position_end, adult_max_cpg_position_end, adult_gene_plus_flanking_cg_methylation_position_end, "methylated_count")
adult_unmethylated_count_cg_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_cpg_position_end, adult_max_cpg_position_end, adult_gene_plus_flanking_cg_methylation_position_end, "unmethylated_count")

adult_percent_chg_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_chg_position, adult_max_chg_position, adult_gene_plus_flanking_chg_methylation_position, "percent")
adult_methylated_chg_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_chg_position, adult_max_chg_position, adult_gene_plus_flanking_chg_methylation_position, "methylated_count")
adult_unmethylated_chg_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_chg_position, adult_max_chg_position, adult_gene_plus_flanking_chg_methylation_position, "unmethylated_count")
adult_percent_chg_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_chg_position_end, adult_max_chg_position_end, adult_gene_plus_flanking_chg_methylation_position_end, "percent")
adult_methylated_count_chg_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_chg_position_end, adult_max_chg_position_end, adult_gene_plus_flanking_chg_methylation_position_end, "methylated_count")
adult_unmethylated_count_chg_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_chg_position_end, adult_max_chg_position_end, adult_gene_plus_flanking_chg_methylation_position_end, "unmethylated_count")

adult_percent_chh_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_chh_position, adult_max_chh_position, adult_gene_plus_flanking_chh_methylation_position, "percent")
adult_methylated_chh_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_chh_position, adult_max_chh_position, adult_gene_plus_flanking_chh_methylation_position, "methylated_count")
adult_unmethylated_chh_methylation_position_tally_dict = position_methylation_stat_tally(adult_min_chh_position, adult_max_chh_position, adult_gene_plus_flanking_chh_methylation_position, "unmethylated_count")
adult_percent_chh_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_chh_position_end, adult_max_chh_position_end, adult_gene_plus_flanking_chh_methylation_position_end, "percent")
adult_methylated_count_chh_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_chh_position_end, adult_max_chh_position_end, adult_gene_plus_flanking_chh_methylation_position_end, "methylated_count")
adult_unmethylated_count_chh_methylation_position_end_tally_dict = position_methylation_stat_tally(adult_min_chh_position_end, adult_max_chh_position_end, adult_gene_plus_flanking_chh_methylation_position_end, "unmethylated_count")

# add length, position, methylated_count, unmethylated_count, median and max data
pupa_cg_position_stats = position_methylation_stats_to_list(pupa_percent_cg_methylation_position_tally_dict, pupa_methylated_count_cg_methylation_position_tally_dict, pupa_unmethylated_count_cg_methylation_position_tally_dict)
pupa_cg_position_end_stats = position_methylation_stats_to_list(pupa_percent_cg_methylation_position_end_tally_dict, pupa_methylated_count_cg_methylation_position_end_tally_dict, pupa_unmethylated_count_cg_methylation_position_end_tally_dict)

pupa_chg_position_stats = position_methylation_stats_to_list(pupa_percent_chg_methylation_position_tally_dict, pupa_methylated_count_chg_methylation_position_tally_dict, pupa_unmethylated_count_chg_methylation_position_tally_dict)
pupa_chg_position_end_stats = position_methylation_stats_to_list(pupa_percent_chg_methylation_position_end_tally_dict, pupa_methylated_count_chg_methylation_position_end_tally_dict, pupa_unmethylated_count_chg_methylation_position_end_tally_dict)

pupa_chh_position_stats = position_methylation_stats_to_list(pupa_percent_chh_methylation_position_tally_dict, pupa_methylated_count_chh_methylation_position_tally_dict, pupa_unmethylated_count_chh_methylation_position_tally_dict)
pupa_chh_position_end_stats = position_methylation_stats_to_list(pupa_percent_chh_methylation_position_end_tally_dict, pupa_methylated_count_chh_methylation_position_end_tally_dict, pupa_unmethylated_count_chh_methylation_position_end_tally_dict)

adult_cg_position_stats = position_methylation_stats_to_list(adult_percent_cg_methylation_position_tally_dict, adult_methylated_cg_methylation_position_tally_dict, adult_unmethylated_cg_methylation_position_tally_dict)
adult_cg_position_end_stats = position_methylation_stats_to_list(adult_percent_cg_methylation_position_end_tally_dict, adult_methylated_count_cg_methylation_position_end_tally_dict, adult_unmethylated_count_cg_methylation_position_end_tally_dict)

adult_chg_position_stats = position_methylation_stats_to_list(adult_percent_chg_methylation_position_tally_dict, adult_methylated_chg_methylation_position_tally_dict, adult_unmethylated_chg_methylation_position_tally_dict)
adult_chg_position_end_stats = position_methylation_stats_to_list(adult_percent_chg_methylation_position_end_tally_dict, adult_methylated_count_chg_methylation_position_end_tally_dict, adult_unmethylated_count_chg_methylation_position_end_tally_dict)

adult_chh_position_stats = position_methylation_stats_to_list(adult_percent_chh_methylation_position_tally_dict, adult_methylated_chh_methylation_position_tally_dict, adult_unmethylated_chh_methylation_position_tally_dict)
adult_chh_position_end_stats = position_methylation_stats_to_list(adult_percent_chh_methylation_position_end_tally_dict, adult_methylated_count_chh_methylation_position_end_tally_dict, adult_unmethylated_count_chh_methylation_position_end_tally_dict)

# plot upstream and gene methylation data in 500 basepair bins (positions -3000 to +3000)
bin_size = 250
start_list = []
end_list = []
i = -3000
while i in range(-3000,3000):
    print(i)
    start_list.append(i)
    end_list.append(i+(bin_size-1))
    i += bin_size

bin_data_list = []
bin_data_label_list = []

for i in range(len(start_list)):
    tmp_bin_data = get_feature_methylation_counts(pupa_gene_plus_flanking_cg_methylation_position, start_list[i], end_list[i])
    bin_data = feature_methylation_counts_to_percent_list(tmp_bin_data)
    bin_data_list.append(bin_data)
    bin_data_label = str(start_list[i])+"."+str(end_list[i])
    bin_data_label_list.append(bin_data_label)

# Create a figure instance
plt.clf()
fig = plt.figure(1, figsize=(9, 6))
# Create an axes instance
ax = fig.add_subplot(111)
ax.set_xticklabels(bin_data_label_list, rotation=40, ha='right')
# Create the boxplot
bp = ax.boxplot(bin_data_list)
# Save the figure
fig.savefig('pupa_upstream_and_gene_methylation_bin_size_250.png', bbox_inches='tight', dpi=300)

# plot number of cpg's analysed between positions -500 and +500
input_position_stats = pupa_chh_position_end_stats
x = []
y = []
z = []
w = []
m = []
u = []
t = []
a = []
for i in range(-5000,5001):
    print(i)
    x.append(i)
    y.append(input_position_stats[1][input_position_stats[0].index(0)+i]) # collect number of values for selected positions
    z.append(input_position_stats[2][input_position_stats[0].index(0)+i]) # collect median values for selected positions
    w.append(input_position_stats[3][input_position_stats[0].index(0)+i])
    m.append(input_position_stats[4][input_position_stats[0].index(0)+i])
    u.append(input_position_stats[5][input_position_stats[0].index(0)+i])
    t.append(input_position_stats[4][input_position_stats[0].index(0)+i]+input_position_stats[5][input_position_stats[0].index(0)+i])
    a.append(input_position_stats[4][input_position_stats[0].index(0)+i]/(input_position_stats[4][input_position_stats[0].index(0)+i]+input_position_stats[5][input_position_stats[0].index(0)+i]))

plt.clf()
plt.plot(x,y, marker='o', linestyle='None', color='b', markersize=2, label='Square')
plt.savefig("pupa_chh_lengthplot_end.png", dpi=300)
plt.clf()
plt.ylim(0, 0.5)
plt.plot(x,z, marker='o', linestyle='None', color='r', markersize=2, label='Square')
plt.savefig("pupa_chh_medianplot_end.png", dpi=300)
plt.clf()
plt.plot(x,w, marker='o', linestyle='None', color='g', markersize=2, label='Square')
plt.savefig("pupa_chh_maxplot_end.png", dpi=300)
plt.clf()
plt.plot(x,t, marker='o', linestyle='None', color='y', markersize=2, label='Square')
plt.savefig("pupa_chh_covplot_end.png", dpi=300)
plt.clf()
plt.plot(x,m, marker='o', linestyle='None', color='m', markersize=2, label='Square')
plt.savefig("pupa_chh_methplot_end.png", dpi=300)
plt.clf()
plt.plot(x,u, marker='o', linestyle='None', color='c', markersize=2, label='Square')
plt.savefig("pupa_chh_unmethplot_end.png", dpi=300)
plt.clf()
plt.ylim(0,0.15)
plt.plot(x,a, marker='o', linestyle='None', color='k', markersize=2, label='Square')
plt.savefig("pupa_chh_methovertotplot_end.png", dpi=300)
plt.clf()
