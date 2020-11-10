
####################################
##                                ##
##  prepareListOfftargetSites.py  ##
##                                ##
####################################

"""
This is the original prepareGeneListsWholeGenome.py file with minor
modifications.
These are meant to allow a choice of parameters through the command line.
They also adjust the assumptions mm10db makes about the structure of the data,
since our use of artificial genomes makes them invalid.
"""

#
# Purpose: identify all offtarget sites in the whole genome
#

#
# Inputs:
#           - one file with all the chromosome sequences, one chromosome per line
#

#
# Output:
#           - one file with all the sites
#


from time import localtime, strftime
import os
import re
import string
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene", type=str, default="hg19",
                        help="The gene to process")
    parser.add_argument("-c", "--chr", type=str, default="chr1",
                        help="The chromosome to process")
    return parser.parse_args()

args = get_args()

# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"


work_dir = os.environ['WORK_DIR']
dir = os.path.join(work_dir, "data/%s/" % args.gene)
dir_seq = os.path.join(dir, "outputs", "mm10db/")

inputFile_chr = "%s.txt" % args.chr
outputFile = "%s_offtargetSites.txt" % args.chr

# This line has been commented out since it is not relevant for artificial
# genomes.
# chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "M", "X", "Y"]


# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq





inFile = open(dir+inputFile_chr,'r')
outFile = open(dir_seq+outputFile,'w')

lineNumber = 0

# For every chromosome
for line_chr in inFile:

    print strftime("%H:%M:%S", localtime())+":\tChromomose number "+ str(lineNumber+1)
    
    # we parse the line and look for forward sequences
    print strftime("%H:%M:%S", localtime())+":\t\tForward-parsing the chromosome."
    match_chr = re.findall(pattern_forward_offsite,line_chr)
    print strftime("%H:%M:%S", localtime())+":\t\tProcessing and saving the parsed sequences."
    if match_chr:
        print "\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites."
        # we save each sequence
        for i in range(0,len(match_chr)):
            outFile.write(match_chr[i][0:20]+"\n")
    else:
        print "\t\t\tWe did not detect any possible off-target sites."

    # we parse the line and look for reverse sequences
    print strftime("%H:%M:%S", localtime())+":\t\tReverse-parsing the chromosome."
    match_chr = re.findall(pattern_reverse_offsite,line_chr)
    print strftime("%H:%M:%S", localtime())+":\t\tProcessing the parsed sequences."
    if match_chr:
        print "\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites."
        # we save each reverse-complement sequence
        for i in range(0,len(match_chr)):
            # we reverse-complement the sequence and count the number of mismatches between the rc and the target
            outFile.write(rc(match_chr[i])[0:20]+"\n")

    lineNumber += 1

inFile.close()
outFile.close()


print "\n"+strftime("%H:%M:%S", localtime())+":\tDone."
