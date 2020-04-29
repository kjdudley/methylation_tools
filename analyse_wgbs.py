####### IMPORT MODULES ########
import matplotlib as mpl
import numpy as np
mpl.use('agg')
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mt_functions import *
from plot_functions import *
#import statistics
import multiprocessing
#import time

### DEFINE PATHS TO DATA FILES ###
path_to_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest"
fasta = "GCF_002156985.1_Harm_1.0_genomic.fa"
gff = "GCF_002156985.1_Harm_1.0_genomic.gff"

#pupa_combined_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/methylationExtractPupaCombine/pupaCombined.CX_report.txt"
#adult_combined_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/200305/methylationExtractAdultCombine/adultCombined.CX_report.txt"
#P4_SSSI_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-SSS1-P4/Jack-methylation-repli-g-SSS1-P4_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
#A7_SSSI_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-SSS1-A7/Jack-methylation-repli-g-SSS1-A7_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
#P4_replig_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-P4/Jack-methylation-repli-g-P4_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
#A7_replig_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-repli-g-A7/Jack-methylation-repli-g-A7_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
#P4_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-P4/Jack-methylation-P4_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
#P5_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-P5/Jack-methylation-P5_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"
#P6_CX_data = "/media/kevlab/projects/helicoverpa_epigenetics/data/20191219/methylation_extract/Jack-methylation-P6/Jack-methylation-P6_R1_merged_trimmed_filtered.deduplicated.CX_report.txt"

#input_file_list = [pupa_combined_CX_data, adult_combined_CX_data, P4_SSSI_CX_data, A7_SSSI_CX_data, P4_replig_CX_data, A7_replig_CX_data, P4_CX_data, P5_CX_data, P6_CX_data]

### IMPORT METHYLATION DATA ###
#pool = multiprocessing.Pool(processes=len(input_file_list))
#sample_methylation_data = list(pool.map(import_methylation_data3, input_file_list))
# save list object for future reuse
import pickle
#with open("pickled_sample_data.txt", "wb") as fp: # pickling
#     pickle.dump(sample_methylation_data, fp)
# uncomment below to reimport
print("*** unpickling! ***")
with open("pickled_sample_data.txt", "rb") as fp: # unpickling
    sample_methylation_data = pickle.load(fp)

contig_list = import_contig_list("/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest/contigs.txt")

for i in contig_list[0:1]:
    print(i)
    for sample in range(len(sample_methylation_data)):
        print(len(sample_methylation_data[sample][i][3]))
        #print("Number of items in dictonary key list:", len(sample_methylation_data[0][i]))
        for j in sample_methylation_data[sample][i][2:5]: # loop through items in dictionary key list, i.e. 6 items
            #print(len(j))
            #print(j[0:2])
           for k in j[0:40]:
                print('{:>3}'.format(k), end="")
           print()

### CREATE A CONTIG SEQUENCE DICTIONARY ###
print("*** creating contig sequence dictionary ***")
print(path_to_data+"/"+fasta)
contig_sequence_dict = import_contig_sequence(contig_list)

### GET GENE ANNOTATION DATA ###
print("*** importing gene annotation data ***")
print(path_to_data+"/"+gff)
gene_annotation = import_gff_data(contig_list, "gene")

### GET NON-GENE COORDINATES ###
print("*** extracting intergene data ***")
intergene_annotation = import_intergene_coordinates(gene_annotation, contig_sequence_dict) # find non genic coordinates in gff format (1-based)

### MARK CPG COORDINATES IN IMPORTED CONTIGS ###
# cpg_coordinates = get_coordinates(contig_sequence_dict, "CG")
# check that it works!
# contig_sequence_dict["NC_014668.1"][cpg_coordinates["NC_014668.1"][0][2]:cpg_coordinates["NC_014668.1"][1][2]]

### EXTRACT STATS FOR GENIC AND INTERGENIC REGIONS ###
# extract stats (i.e. feature size, cpg count, cpg percent, cpg coordinates)
# contig_gene_stats_dict = get_feature_stats(cpg_coordinates, gene_annotation, contig_sequence_dict)
# contig_intergene_stats_dict = get_feature_stats(cpg_coordinates, intergene_annotation, contig_sequence_dict)

### ANALYSE SAMPLES ###
required_coverage = 20
cg_context = "CG"
chg_context = "CHG"
chh_context = "CHH"
print("*** starting sample processing ***")
sample_methylation_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_methylation_dict[str(sample_number)+"_"+cg_context+"_contig_methylation_stats_dict"] = get_sample_methylation_stats(sample_methylation_data[sample], required_coverage, cg_context)
    sample_methylation_dict[str(sample_number)+"_"+chg_context+"_contig_methylation_stats_dict"] = get_sample_methylation_stats(sample_methylation_data[sample], required_coverage, chg_context)
    sample_methylation_dict[str(sample_number)+"_"+chh_context+"_contig_methylation_stats_dict"] = get_sample_methylation_stats(sample_methylation_data[sample], required_coverage, chh_context)

# group data by gene and intergene (i.e. feature) coordinates
print("*** grouping data by feature ***")
sample_feature_methylation_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_methylation_stats_dict"] = contig_to_feature(sample_methylation_dict[str(sample_number)+"_"+cg_context+"_contig_methylation_stats_dict"], gene_annotation)
    sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_methylation_stats_dict"] = contig_to_feature(sample_methylation_dict[str(sample_number)+"_"+chg_context+"_contig_methylation_stats_dict"], gene_annotation)
    sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_methylation_stats_dict"] = contig_to_feature(sample_methylation_dict[str(sample_number)+"_"+chh_context+"_contig_methylation_stats_dict"], gene_annotation)
    sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_intergene_methylation_stats_dict"] = contig_to_feature(sample_methylation_dict[str(sample_number)+"_"+cg_context+"_contig_methylation_stats_dict"], intergene_annotation)
    sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_intergene_methylation_stats_dict"] = contig_to_feature(sample_methylation_dict[str(sample_number)+"_"+chg_context+"_contig_methylation_stats_dict"], intergene_annotation)
    sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_intergene_methylation_stats_dict"] = contig_to_feature(sample_methylation_dict[str(sample_number)+"_"+chh_context+"_contig_methylation_stats_dict"], intergene_annotation)

# add flanking regions to gene coordinate data (upstream and downstream)
print("*** adding flanking regions ***")
sample_gene_plus_flanking_methylation_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_stats_dict"] = \
add_upstream_downstream_regions(sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_methylation_stats_dict"], sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_intergene_methylation_stats_dict"])
    sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_stats_dict"] = \
add_upstream_downstream_regions(sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_methylation_stats_dict"], sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_intergene_methylation_stats_dict"])
    sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_stats_dict"] = \
add_upstream_downstream_regions(sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_methylation_stats_dict"], sample_feature_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_intergene_methylation_stats_dict"])

# convert coordinates to positions
print("*** converting coordinates to position values ***")
sample_gene_plus_flanking_methylation_position_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"] = \
gene_coordinate_to_position(sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_stats_dict"], "start")
    sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"] = \
gene_coordinate_to_position(sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_stats_dict"], "start")
    sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_start_stats_dict"] = \
gene_coordinate_to_position(sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_stats_dict"], "start")
    sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"] = \
gene_coordinate_to_position(sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_stats_dict"], "end")
    sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"] = \
gene_coordinate_to_position(sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_stats_dict"], "end")
    sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_end_stats_dict"] = \
gene_coordinate_to_position(sample_gene_plus_flanking_methylation_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_stats_dict"], "end")

# find max and min positions
print("*** finding max and min values ***")
sample_gene_plus_flanking_methylation_max_position_stats_data_dict = {}
sample_gene_plus_flanking_methylation_min_position_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"] = \
get_max_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"])
    sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"] = \
get_max_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"])
    sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"] = \
get_max_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_start_stats_dict"])
    sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"] = \
get_min_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"])
    sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"] = \
get_min_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"])
    sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"] = \
get_min_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"])
    sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"] = \
get_max_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"])
    sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"] = \
get_max_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"])
    sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"] = \
get_max_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_end_stats_dict"])
    sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"] = \
get_min_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"])
    sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"] = \
get_min_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"])
    sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"] = \
get_min_cpg_position(sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"])

# create percent cpg dictionary entry for each position between min and max
print("*** counting cpg methylation position calls ***")
sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_percent_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "percent")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_methylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "methylated_count")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_unmethylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "unmethylated_count")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_percent_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "percent")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_methylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "methylated_count")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_unmethylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "unmethylated_count")

print("*** counting chg methylation position calls ***")
sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_percent_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "percent")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_methylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "methylated_count")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_unmethylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "unmethylated_count")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_percent_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "percent")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_methylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "methylated_count")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_unmethylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "unmethylated_count")

print("*** counting chh methylation position calls ***")
sample_gene_plus_flanking_cpg_methylation_max_to_min_position_stats_data_dict = {}
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_start_percent_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "percent")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_start_methylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "methylated_count")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_start_unmethylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_start_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_start_stats_dict"], "unmethylated_count")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_end_percent_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "percent")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_end_methylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "methylated_count")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_end_unmethylated_count_stats_data_dict"] = position_methylation_stat_tally(\
sample_gene_plus_flanking_methylation_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_min_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_max_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_position_end_stats_dict"], \
sample_gene_plus_flanking_methylation_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_position_end_stats_dict"], "unmethylated_count")

# add length, position, methylated_count, unmethylated_count, median and max data
print("*** calculating cpg position methylation statistics ***")
sample_gene_plus_flanking_cpg_methylation_max_to_min_position_summary_stats_data_list = []
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_summary_stats_data_list.append(position_methylation_stats_to_list(\
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_percent_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_methylated_count_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_unmethylated_count_stats_data_dict"]))
    sample_gene_plus_flanking_cpg_methylation_max_to_min_position_summary_stats_data_list.append(position_methylation_stats_to_list(\
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_percent_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_methylated_count_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+cg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_unmethylated_count_stats_data_dict"]))

print("*** calculating chg position methylation statistics ***")
sample_gene_plus_flanking_chg_methylation_max_to_min_position_summary_stats_data_list = []
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_summary_stats_data_list.append(position_methylation_stats_to_list(\
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_percent_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_methylated_count_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_start_unmethylated_count_stats_data_dict"]))
    sample_gene_plus_flanking_chg_methylation_max_to_min_position_summary_stats_data_list.append(position_methylation_stats_to_list(\
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_percent_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_methylated_count_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chg_context+"_gene_plus_flanking_methylation_max_to_min_position_end_unmethylated_count_stats_data_dict"]))

print("*** calculating chh position methylation statistics ***")
sample_gene_plus_flanking_chh_methylation_max_to_min_position_summary_stats_data_list = []
for sample in range(len(sample_methylation_data)):
    sample_number = sample+1
    print("*** Processing sample", sample_number, "***")
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_summary_stats_data_list.append(position_methylation_stats_to_list(\
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_start_percent_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_start_methylated_count_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_start_unmethylated_count_stats_data_dict"]))
    sample_gene_plus_flanking_chh_methylation_max_to_min_position_summary_stats_data_list.append(position_methylation_stats_to_list(\
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_end_percent_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_end_methylated_count_stats_data_dict"], \
sample_gene_plus_flanking_methylation_max_to_min_position_stats_data_dict[str(sample_number)+"_"+chh_context+"_gene_plus_flanking_methylation_max_to_min_position_end_unmethylated_count_stats_data_dict"]))

# plot upstream and gene methylation data in 500 basepair bins (positions -3000 to +3000)
#bin_size = 250
#start_list = []
#end_list = []
#i = -3000
#while i in range(-3000,3000):
#    print(i)
#    start_list.append(i)
#    end_list.append(i+(bin_size-1))
#    i += bin_size

#bin_data_list = []
#bin_data_label_list = []

#for i in range(len(start_list)):
#    tmp_bin_data = get_feature_methylation_counts(pupa_gene_plus_flanking_cg_methylation_position, start_list[i], end_list[i])
#    bin_data = feature_methylation_counts_to_percent_list(tmp_bin_data)
#    bin_data_list.append(bin_data)
#    bin_data_label = str(start_list[i])+"."+str(end_list[i])
#    bin_data_label_list.append(bin_data_label)

# Create a figure instance
#plt.clf()
#fig = plt.figure(1, figsize=(9, 6))
# Create an axes instance
#ax = fig.add_subplot(111)
#ax.set_xticklabels(bin_data_label_list, rotation=40, ha='right')
# Create the boxplot
#bp = ax.boxplot(bin_data_list)
# Save the figure
#fig.savefig('pupa_upstream_and_gene_methylation_bin_size_250.png', bbox_inches='tight', dpi=300)
# plot number of cpg's analysed between positions -500 and +500

#all_sample_position_stats_list = [pupa_cg_position_stats, pupa_cg_position_end_stats, pupa_chg_position_stats, pupa_chg_position_end_stats, pupa_chh_position_stats, pupa_chh_position_end_stats, adult_cg_position_stats, adult_cg_position_end_stats, adult_chg_position_stats, adult_chg_position_end_stats, adult_chh_position_stats, adult_chh_position_end_stats, P4_SSSI_cg_position_stats, P4_SSSI_cg_position_end_stats, P4_SSSI_chg_position_stats, P4_SSSI_chg_position_stats, P4_SSSI_chh_position_stats, P4_SSSI_chh_position_end_stats, A7_SSSI_cg_position_stats, A7_SSSI_cg_position_end_stats, A7_SSSI_chg_position_stats, A7_SSSI_chg_position_end_stats, A7_SSSI_chh_position_stats, A7_SSSI_chh_position_end_stats, P4_replig_cg_position_stats, P4_replig_cg_position_end_stats, P4_replig_chg_position_stats, P4_replig_chg_position_end_stats, P4_replig_chh_position_stats, P4_replig_chh_position_end_stats, A7_replig_cg_position_stats, A7_replig_cg_position_end_stats, A7_replig_chg_position_stats, A7_replig_chg_position_end_stats, A7_replig_chh_position_stats, A7_replig_chh_position_end_stats]

### PLOT DATA ###
#print("*** plotting cpg data ***")
#all_cpg_position_stats_plot_data = []
#sample_count = 0
#for i in sample_gene_plus_flanking_cpg_methylation_max_to_min_position_summary_stats_data_list:    
#    print("*** processing sample", i, "***")
#    sample_position_stats_plot_data = []
#    x = []
#    y = []
#    z = []
#    w = []
#    m = []
#    u = []
#    t = []
#    a = []
#    for j in range(-3000,3001):
#        x.append(j)
#        y.append(i[1][i[0].index(0)+j]) # collect number of values for selected positions
#        z.append(i[2][i[0].index(0)+j]) # collect median values for selected positions
#        w.append(i[3][i[0].index(0)+j])
#        m.append(i[4][i[0].index(0)+j])
#        u.append(i[5][i[0].index(0)+j])
#        t.append(i[4][i[0].index(0)+j]+i[5][i[0].index(0)+j])
#        a.append(i[4][i[0].index(0)+j]/(i[4][i[0].index(0)+j]+i[5][i[0].index(0)+j]))
#    sample_cpg_position_stats_plot_data = [x, y, z, w, m, u, t, a]
#    all_cpg_position_stats_plot_data.append(sample_position_stats_plot_data)

sample_count = 0
for i in range(len(all_cpg_position_stats_plot_data[0])):
    sample_count += 1
    print("*** plotting percent cpg methylation data for sample", sample_count, "***")
    x_axis = all_cpg_position_stats_plot_data[i][0]
    y_axis = all_cpg_position_stats_plot_data[i][7]
    #for k in range(len(x_axis)):
    #    print(x_axis[k], y_axis[k], sep=":", end=" ")
    print("number of datapoints:", len(y_axis))
    print("min percent:", min(y_axis)*100)
    print("max percent:", max(y_axis)*100)
    print("mean percent:", sum(y_axis)/len(y_axis))
    plt.clf()
    plt.plot(x_axis, y_axis, marker='o', linestyle='None', color='k', markersize=2, label='Square')
    plt.savefig(f'{str(i+1)}_percent_cpg_methylation.png')
    plt.clf()

#plt.clf()
#plt.plot(x,y, marker='o', linestyle='None', color='b', markersize=2, label='Square')
#plt.savefig("pupa_chh_lengthplot_end.png", dpi=300)
#plt.clf()
#plt.ylim(0, 0.5)
#plt.plot(x,z, marker='o', linestyle='None', color='r', markersize=2, label='Square')
#plt.savefig("pupa_chh_medianplot_end.png", dpi=300)
#plt.clf()
#plt.plot(x,w, marker='o', linestyle='None', color='g', markersize=2, label='Square')
#plt.savefig("pupa_chh_maxplot_end.png", dpi=300)
#plt.clf()
#plt.plot(x,t, marker='o', linestyle='None', color='y', markersize=2, label='Square')
#plt.savefig("pupa_chh_covplot_end.png", dpi=300)
#plt.clf()
#plt.plot(x,m, marker='o', linestyle='None', color='m', markersize=2, label='Square')
#plt.savefig("pupa_chh_methplot_end.png", dpi=300)
#plt.clf()
#plt.plot(x,u, marker='o', linestyle='None', color='c', markersize=2, label='Square')
#plt.savefig("pupa_chh_unmethplot_end.png", dpi=300)
#plt.clf()
#plt.ylim(0,0.15)
#plt.plot(x,a, marker='o', linestyle='None', color='k', markersize=2, label='Square')
#plt.savefig("pupa_chh_methovertotplot_end.png", dpi=300)
#plt.clf()
