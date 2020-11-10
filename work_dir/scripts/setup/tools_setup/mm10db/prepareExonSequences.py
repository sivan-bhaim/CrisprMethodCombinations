
##################################
##                              ##
##   prepareExonSequences.py    ##
##                              ##
##################################

"""
This is the original prepareExonSequences.py file with minor modifications.
These are meant to allow a choice of parameters through the command line.
They also adjust the assumptions mm10db makes about the structure of the data,
since our use of artificial genomes makes them invalid.
"""

#
# Purpose: preparing the sequences for all exons of interest
#

#
# Inputs:
#           - one file with all the chromosome sequences, one chromosome per line
#           - the list of exons of interest
#

#
# Output:
#           - one file with all the exon sequences, one exon per line
#           - note: we include 17bp before and after the exon
#


from time import localtime, strftime, sleep
import re
import ast
import sys
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene", type=str, default="hg19",
                        help="The gene to process")
    parser.add_argument("-c", "--chr", type=str, default="chr1",
                        help="The chromosome to process")
    parser.add_argument("-l", "--list", type=str, default="",
                        help="The list file")
    return parser.parse_args()

args = get_args()

work_dir = os.environ['WORK_DIR']
dir = os.path.join(work_dir, "data/%s/" % args.gene)
dir_seq = os.path.join(dir, "outputs", "mm10db/")
dir_list = dir_seq
chr_ = "%s.txt" % args.chr
list_ = "%s_exonList.txt" % args.chr
exons_ = "%s_exon_sequences.txt" % args.chr

if args.list == "":
    print "Extracting exon sequences for the whole genome"

else:
    list_ = list_[:-4]+"_"+sys.argv[1]+".txt"
    print "Extracting sequence for exons listed in "+list_
    if os.path.isfile(dir_list+list_)==False:
        print "Error: file "+list_+" does not exist."
        quit()
    exons_ = exons_[:-4]+"_"+sys.argv[1]+".txt"


padding = 17 # digestion site must be on exon, but the whole target does not have to

# This line has been commented out since it is not relevant for artificial
# genomes.
# chr_offset = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7, "9": 8, "10": 9, "11": 10, "12": 11, "13": 12, "14": 13, "15": 14, "16": 15, "17": 16, "18": 17, "19": 18, "M": 19, "X": 20, "Y": 21}

pattern = r".+\tchr(.+)\t[\+-]\t(\d+)\t(\d+)"

# First, we load all the chromosome sequences
print strftime("%H:%M:%S", localtime())+": Loading chromosome sequences from "+chr_
inFile = open(dir+chr_,'r')
chrSeq = inFile.readlines()
inFile.close()


# Then, we load all the list of exons
print strftime("%H:%M:%S", localtime())+": Loading the exon list from "+list_
inFile = open(dir_list+list_,'r')
exonList = inFile.readlines()
inFile.close()


# Finally, we can extract (and save) the sequence for each exon (with padding on each side)
print strftime("%H:%M:%S", localtime())+": Extracting and saving the exon sequences"
outFile = open(dir_seq+exons_,'w')
for i in range(0,len(exonList)):
    match = re.search(pattern,exonList[i])
    if match:
        # We have hardcoded chr=0 since we have no offsets in the artificial
        # genome.
        # chr = chr_offset[match.group(1)]
        chr = 0
        start = ast.literal_eval(match.group(2))-1
        end = ast.literal_eval(match.group(3))
        outFile.write(chrSeq[chr][start-padding:end+padding].upper()+"\n")
    else:
        print "Problem? "+exonList[i].rstrip()


print strftime("%H:%M:%S", localtime())+": Done."



#######################
##                   ##
##    End of File    ##
##                   ##
#######################
