
##################################
##                              ##
##      createListExons.py      ##
##                              ##
##################################

"""
This is the original createListExons.py file with minor modifications.
These are meant to allow a choice of parameters through the command line.
"""

#
# Purpose: identifying all exons that may contain targets of interest
# If a gene has two or more isoforms, this means only common exons.
# Otherwise, this means all exons.
#

#
# Inputs:
#           - refGene file, that contains all genes
#           - optionally, a list of genes for which we want to extract exons
#

#
# Output:
#           - an "exons of interest" file, where each line is an exon
#


import re
import sys
import os
import ast
import argparse
from time import localtime, strftime, sleep
from subprocess import call


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene", type=str, default="xu",
                        help="The gene to process")
    parser.add_argument("-c", "--chr", type=str, default="chr1",
                        help="The chromosome to process")
    parser.add_argument("-l", "--list", type=str, default="",
                        help="The list file")
    return parser.parse_args()

args = get_args()


work_dir = os.environ['WORK_DIR']
dir = os.path.join(work_dir, "data/%s/" % args.gene)
mm10_dir = os.path.join(dir, "outputs", "mm10db")
in_ = "%s_refGene.txt" % args.chr
out_ = "%s_exonList.txt" % args.chr

pattern = r"\d+\t(N[MR]_\d+)\t(chr[MXY\d]+)\t([\+-])\t\d+\t\d+\t(\d+)\t(\d+)\t\d+\t([\d,]+)\t([\d,]+)\t0\t(.+)\t(cmpl|unk|incmpl)\t(cmpl|unk|incmpl)"

if args.list == "":
    print "Creating a list of exons for the whole genome"

else:
    filename = args.list+".txt"
    print "Creating a list of exons for genes listed in "+filename
    if os.path.isfile(dir+filename)==False:
        print "Error: file "+filename+" does not exist."
        quit()
    genes = r""
    inFile = open(os.path.join(mm10_dir, filename),'r')
    for line in inFile:
        genes += line.rstrip()+"|"

    pattern = r"\d+\t(N[MR]_\d+)\t(chr[MXY\d]+)\t([\+-])\t\d+\t\d+\t(\d+)\t(\d+)\t\d+\t([\d,]+)\t([\d,]+)\t0\t("
    pattern += genes[:-1]
    pattern += r")\t(cmpl|unk|incmpl)\t(cmpl|unk|incmpl)"

    out_ = out_[:-4]+"_"+filename



print strftime("%H:%M:%S", localtime())+": Loading data from "+in_


savedData = []
nb = 0
genes = dict()

# We start by reading the refGene file
inFile = open(os.path.join(dir, in_),'r')
for line in inFile:
    match = re.search(pattern,line.rstrip())
    if match:
        # For each gene isoform, we save the ID, the chromosome and direction, the exon positions, and the gene name
        savedData.append(match.group(1)+"\t"+match.group(2)+"\t"+match.group(3)+"\t"+match.group(6)+"\t"+match.group(7)+"\t"+match.group(8)+"\t"+match.group(4)+"\t"+match.group(5))
        currentGene = match.group(8)
        
        # We also save the position of this data into a dictionnary. This will help detecting genes with multiple isoforms
        if currentGene in genes:
            genes[currentGene].append(nb)
        else:
            genes[currentGene] = []
            genes[currentGene].append(nb)
        nb+=1

# Next, we process all this data
print strftime("%H:%M:%S", localtime())+": Processing data, and saving to "+out_

out_path = os.path.join(mm10_dir, out_)
outFile = open(out_path,'w')

for currentGene in genes:
    nbIsoforms = len(genes[currentGene])

    # If the gene has a unique isoform, all exons are interesting
    if nbIsoforms==1:
        # We locate the position of the saved data
        offset = genes[currentGene][0]
        # We extract the data
        split_data = savedData[offset].split("\t")
        chr = split_data[1]
        direction = split_data[2]
        name = split_data[5]
        exon_starts = split_data[3].rstrip(",").split(",")
        exon_ends = split_data[4].rstrip(",").split(",")
        cds_start = split_data[6]
        cds_end = split_data[7]
        # We save the data, exon by exon
        for i in range(0,len(exon_starts)):
            exon_data = name+"\t"+chr+"\t"+direction+"\t"+exon_starts[i]+"\t"+exon_ends[i]+"\t"+cds_start+"\t"+cds_end+"\tNo_Isoform"
            outFile.write(exon_data+"\n")

    # If the gene has two or more isoforms, only the common exons are interesting
    else:
        exon_overlap = dict()
        # First, we prepare the CDS information
        cds_start = ""
        cds_end = ""
        cds_status = ""
        for j in range(0,nbIsoforms):
            offset = genes[currentGene][j]
            # We extract the data
            split_data = savedData[offset].split("\t")
            # For the first isoform, we simply initialise
            if j==0:
                cds_start = split_data[6]
                cds_end = split_data[7]
                cds_status = "Unique_CDS"
            # For the next ones, we compare the CDS information
            else:
                # If the CDS is the same, nothing happens
                if cds_start == split_data[6] and cds_end == split_data[7]:
                    continue
                # Otherwise, we save the intersection of the two regions
                else:
                    old_CDS_start = ast.literal_eval(cds_start)
                    old_CDS_end = ast.literal_eval(cds_end)
                    new_CDS_start = ast.literal_eval(split_data[6])
                    new_CDS_end = ast.literal_eval(split_data[7])
                    cds_start = str(max(old_CDS_start,new_CDS_start))
                    cds_end = str(min(old_CDS_end,new_CDS_end))
                    cds_status = "CDS_Intersection"
        # We locate the position of the saved data for each isoform
        for j in range(0,nbIsoforms):
            offset = genes[currentGene][j]
            # We extract the data
            split_data = savedData[offset].split("\t")
            chr = split_data[1]
            direction = split_data[2]
            name = split_data[5]
            exon_starts = split_data[3].rstrip(",").split(",")
            exon_ends = split_data[4].rstrip(",").split(",")
            # We add the data, exon by exon, to the dictionnary
            for i in range(0,len(exon_starts)):
                exon_data = name+"\t"+chr+"\t"+direction+"\t"+exon_starts[i]+"\t"+exon_ends[i]+"\t"+cds_start+"\t"+cds_end+"\t"+cds_status
                if exon_data in exon_overlap:
                    exon_overlap[exon_data] += 1
                else:
                    exon_overlap[exon_data] = 1
        # Now, only exons that are present all the time are saved
        # Note: we only accept total overlap
        # If, in one isoform, the start or end codon is different, we reject all versions of this exon
        for ex in exon_overlap:
            if exon_overlap[ex] == nbIsoforms:
                outFile.write(ex+"\n")

outFile.close()

temp_path = os.path.join(mm10_dir, "temp.txt")
cmd = "sort -k 2,2 -k4,4n "+out_path+" > "+temp_path+" && mv "+temp_path+" "+out_path
call([cmd],shell=True)

print strftime("%H:%M:%S", localtime())+": Done."


#######################
##                   ##
##    End of File    ##
##                   ##
#######################
