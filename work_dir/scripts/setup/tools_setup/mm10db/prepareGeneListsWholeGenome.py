
#######################################
##                                   ##
##  prepareGeneListsWholeGenome.py   ##
##                                   ##
#######################################

"""
This is the original prepareListOfftargetSites.py file with minor modifications.
These are meant to allow a choice of parameters through the command line.
"""

#
# Purpose: preparing lists of genes covering the whole genome
#

#
# Output:
#           - several lists, with at most 20 genes in each
#


import argparse
from time import localtime, strftime
import re
import ast

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene", type=str, default="hg19",
                        help="The gene to process")
    parser.add_argument("-c", "--chr", type=str, default="chr1",
                        help="The chromosome to process")
    return parser.parse_args()

args = get_args()

work_dir = os.environ['WORK_DIR']
dir = os.path.join(work_dir, "data/%s/" % args.gene)
mm10_dir = dir + "mm10db/"
refGene = "refGene.txt"

pattern = r"\d+\t(N[MR]_\d+)\t(chr[MXY\d]+)\t([\+-])\t\d+\t\d+\t(\d+)\t(\d+)\t\d+\t([\d,]+)\t([\d,]+)\t0\t(.+)\t(cmpl|unk|incmpl)\t(cmpl|unk|incmpl)"

max_number_of_genes = 25


print strftime("%H:%M:%S", localtime())+": Preparing lists covering the whole genome."


# load all genes from reference genome

print strftime("%H:%M:%S", localtime())+":\tLoading genes from reference genome."

refSeqGenes = []
inFile = open(dir+refGene,'r')
for line in inFile:
    match = re.search(pattern,line.rstrip())
    if match:
        geneName = match.group(8)
        if geneName not in refSeqGenes:
            refSeqGenes.append(geneName)
inFile.close()

print strftime("%H:%M:%S", localtime())+":\t"+str(len(refSeqGenes))+" are yet to be listed."



# Finally, we create the lists

print strftime("%H:%M:%S", localtime())+":\tSaving new lists."

list_index = 0
for i in range(0,len(refSeqGenes)):
    if(i%max_number_of_genes==0):
        list_index += 1
        filename = "list_%03d.txt" % list_index
        if i>0:
            outFile.close()
        outFile = open(mm10_dir+filename,'w')
    outFile.write(refSeqGenes[i]+"\n")
outFile.close()


print strftime("%H:%M:%S", localtime())+": Done. "+str(list_index)+" lists created."



#######################
##                   ##
##    End of File    ##
##                   ##
#######################
