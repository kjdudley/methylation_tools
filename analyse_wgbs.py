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
cx_report = "adultCombined.CX_report.txt"
sample_CX_data_file = path_to_data+"/"+cx_report
fasta = "GCF_002156985.1_Harm_1.0_genomic.fa"
gff = "GCF_002156985.1_Harm_1.0_genomic.gff"

### CREATE A CONTIG SEQUENCE DICTIONARY ###

#contig_list = ["NC_014668.1","NW_018396379.1","NW_018396386.1","NW_018396377.1","NW_018396384.1"]
contig_list = import_contig_list("/media/kevlab/projects/helicoverpa_epigenetics/data/200305/coverageStatsTest/contigs_test.txt")
contig_sequence_dict = import_contig_sequence(contig_list)

print(len(contig_sequence_dict), "contigs to be analysed!")

for i in contig_sequence_dict:
    print(i,"\t", contig_sequence_dict[i][0:101]) # print first 100 nucleotides of each contig

### GET GENE ANNOTATION DATA ###

gene_annotation = import_gff_data(contig_list, "gene")

### GET NON-GENE COORDINATES ###

intergene_annotation = import_intergene_coordinates(gene_annotation, contig_sequence_dict) # find non genic coordinates in gff format (1-based)

### IMPORT METHYLATION DATA ###

adult_combined_methylation_dict = import_methylation_data(contig_list, sample_CX_data_file)

### MARK CPG COORDINATES IN IMPORTED CONTIGS ###

cpg_coordinates = get_coordinates(contig_sequence_dict, "CG")
# check that it works!
# contig_sequence_dict["NC_014668.1"][cpg_coordinates["NC_014668.1"][0][2]:cpg_coordinates["NC_014668.1"][1][2]]

### EXTRACT SPECIFIC METRICS ###

# extract gene stats (i.e. size, cpg count, cpg percent, cpg coordinates)

contig_gene_cpg_count_dict = {}
contig_gene_cpg_coords_dict = {}
contig_gene_size_dict = {}
contig_gene_cpg_percent_dict = {}

for i in contig_sequence_dict:
    print(i)
    number_of_cpgs = []
    cpg_coords = []
    gene_cpg_coords_dict = {}
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
        gene_cpg_coords_dict[gene_start_coordinate] = gene_cpg_coords
        size_of_gene.append(gene_end_coordinate-gene_start_coordinate)
        gene_percent_cpg.append(100*(gene_cpg_count/(gene_end_coordinate-gene_start_coordinate)))
        print("gene_cpg_count:", gene_cpg_count)
    contig_gene_cpg_count_dict[i] = number_of_cpgs
    contig_gene_cpg_coords_dict[i] = gene_cpg_coords_dict
    contig_gene_size_dict[i] = size_of_gene
    contig_gene_cpg_percent_dict[i] = gene_percent_cpg

### ANALYSE SAMPLES ###

contig_cpg_coverage_threshold_dict = {}
contig_cpg_coordinate_dict = {}
contig_cpg_percent_methylation_dict = {}
contig_cpg_counts_methylation_dict = {}
required_coverage = 20

for i in contig_gene_cpg_coords_dict: # contig name
    print(i)
    contig_coordinate = []
    contig_percent_data = []
    contig_cumulative_methylated = 0
    contig_cumulative_unmethylated = 0
    contig_pass_coverage_threshold = 0
    contig_fail_coverage_threshold = 0
    for j in contig_gene_cpg_coords_dict[i]: # gene start coordinate
        for k in contig_gene_cpg_coords_dict[i][j]: # cpg coordinate
            #coordinate = adult_combined_methylation_dict[i][0][adult_combined_methylation_dict[i][0].index(k)] # not required as coordinate equals k!
            coordinate = k
            methylated = adult_combined_methylation_dict[i][2][adult_combined_methylation_dict[i][0].index(k)]
            unmethylated = adult_combined_methylation_dict[i][3][adult_combined_methylation_dict[i][0].index(k)]
            total = methylated+unmethylated
            context = adult_combined_methylation_dict[i][4][adult_combined_methylation_dict[i][0].index(k)]
            if total >= required_coverage:
                contig_coordinate.append(coordinate)
                percent = 100*(methylated/total)
                contig_percent_data.append(percent)
                contig_cumulative_methylated += adult_combined_methylation_dict[i][2][adult_combined_methylation_dict[i][0].index(k)]
                contig_cumulative_unmethylated += adult_combined_methylation_dict[i][3][adult_combined_methylation_dict[i][0].index(k)]
                contig_pass_coverage_threshold += 1
            else:
                contig_fail_coverage_threshold += 1
    contig_cpg_coverage_threshold_dict[i] = [contig_pass_coverage_threshold, contig_fail_coverage_threshold]
    contig_cpg_coordinate_dict[i] = contig_coordinate
    contig_cpg_percent_methylation_dict[i] = contig_percent_data
    contig_totals_data = [contig_cumulative_methylated, contig_cumulative_unmethylated]
    contig_cpg_counts_methylation_dict[i] = contig_totals_data
    if contig_cumulative_methylated+contig_cumulative_unmethylated > 0: # discounts contigs without any CpG's that pass the threshold!
        print(contig_cpg_percent_methylation_dict[i])
        print(100*(pass_coverage_threshold/(pass_coverage_threshold+fail_coverage_threshold)), "percent of CpG dinucleotides pass coverage threshold!")
        print("Total methylated CpG's:  ", cumulative_methylated)
        print("Total unmethylated CpG's:", cumulative_unmethylated)
        print("Overall CpG methylation for passing CpG dinucleotides is", 100*(cumulative_methylated/(cumulative_methylated+cumulative_unmethylated)))
    else:
        print("No data for this contig!")

print(contig_counts_methylation_dict)
#print(contig_cpg_percent_methylation_dict)

contig_gene_cpg_position_methylation_dict = gene_coordinate_to_position(contig_gene_cpg_coordinate_methylation_dict)

contig_max_cpg_position = {}
for i in contig_gene_cpg_position_methylation_dict:
    max_cpg_position = 0
    for j in contig_gene_cpg_position_methylation_dict[i]:
        if len(contig_gene_cpg_position_methylation_dict[i][j][3]) > 0:
            if max(contig_gene_cpg_position_methylation_dict[i][j][3]) > max_cpg_position:
                max_cpg_position = max(contig_gene_cpg_position_methylation_dict[i][j][3])
    contig_max_cpg_position[i] = max_cpg_position

contig_min_cpg_position = {}
for i in contig_gene_cpg_position_methylation_dict:
    min_cpg_position = contig_max_cpg_position[i]
    for j in contig_gene_cpg_position_methylation_dict[i]:
        if len(contig_gene_cpg_position_methylation_dict[i][j][3]) > 0:
            if min(contig_gene_cpg_position_methylation_dict[i][j][3]) < min_cpg_position:
                min_cpg_position = min(contig_gene_cpg_position_methylation_dict[i][j][3])
    contig_min_cpg_position[i] = min_cpg_position

position_percent_cpg_list_dict = {}
for i in range(0, max(contig_max_cpg_position.values())+1):
    position_percent_cpg_list_dict[i] = []

for i in contig_gene_cpg_position_methylation_dict:
    for j in contig_gene_cpg_position_methylation_dict[i]:
        if len(contig_gene_cpg_position_methylation_dict[i][j][3]) > 0:
            #print("TRUE")
            for k in range(len(contig_gene_cpg_position_methylation_dict[i][j][3])):
                #print(len(contig_gene_cpg_position_methylation_dict[i][j][3]))
                print("position:", contig_gene_cpg_position_methylation_dict[i][j][3][k], "percent:", contig_gene_cpg_position_methylation_dict[i][j][4][k])
                position_list = position_percent_cpg_list_dict[contig_gene_cpg_position_methylation_dict[i][j][3][k]]
                position_list.append(contig_gene_cpg_position_methylation_dict[i][j][4][k])
                position_percent_cpg_list_dict[contig_gene_cpg_position_methylation_dict[i][j][3][k]] = position_list

position_length = []
for i in position_percent_cpg_list_dict:
     position_length.append(len(position_percent_cpg_list_dict[i]))

position_median_percent = []
for i in position_percent_cpg_list_dict:
   if len(position_percent_cpg_list_dict[i]) > 0:
      position_median_percent.append(statistics.median(position_percent_cpg_list_dict[i]))
   else:
      position_median_percent.append(-1)

### PLOT DATA ###

## Boxplot
data_to_plot = []
for i in contig_cpg_percent_methylation_dict:
    name = tuple(contig_cpg_percent_methylation_dict[i])
    data_to_plot.append(name)

# Create a figure instance
fig = plt.figure(1, figsize=(9, 6))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(data_to_plot)

# Save the figure
fig.savefig('fig2.png', bbox_inches='tight')

## Histogram
hist_data = []
for i in contig_cpg_percent_methylation_dict:
    name = list(contig_cpg_percent_methylation_dict[i])
    hist_data.append(name)

# create bins
bin = []
for i in range(0,101):
    bin.append(i)

# assign bin values
bin_values = []

for i in bin:
    for j in hist_data[1]:
        if i == int(round(j, 0)):
            binned_values.append(i)

# create bin frequency dict
counted = count_elements(binned_values)

# create histogram
ascii_histogram(binned_values)

## Barplots
plt.clf()
x = np.arange(len(position_length[0:5000]))
y = position_length[0:5000]
plt.bar(x,y)
plt.savefig("lengthbarplot.png")
plt.clf()

y = position_median_percent[0:5000]
plt.bar(x,y)
plt.savefig("medianbarplot.png")
plt.clf()

