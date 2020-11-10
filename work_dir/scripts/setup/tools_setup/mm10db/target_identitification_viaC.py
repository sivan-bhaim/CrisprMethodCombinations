
##################################
##                              ##
##  target_identitification.py  ##
##                              ##
##################################

"""
This is the original target_identitification_viaC.py file with minor
modifications.
These are meant to allow a choice of parameters through the command line.
They also adjust the assumptions mm10db makes about the structure of the data,
since our use of artificial genomes makes them invalid.
"""

#
# Purpose: identify all potential targets in exons
#

#
# Inputs:
#           - one file with all the exon sequences, one exon per line
#           - the list of exons of interest
#

#
# Output:
#           - selected targets
#


from subprocess import call
import re
import ast
from time import localtime, strftime
from math import log
import glob
import os
import string
import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene", type=str, default="hg19",
                        help="The gene to process")
    parser.add_argument("-c", "--chr", type=str, default="chr1",
                        help="The chromosome to process")
    parser.add_argument("-l", "--list", type=str, default="all",
                        help="The genes list file")
    parser.add_argument("--tc", type=int, default=1,
                        help="nb_threads_C")
    parser.add_argument("--tb", type=int, default=1,
                        help="nb_threads_Bowtie")
    return parser.parse_args()

args = get_args()


# Defining the patterns used to detect sequences
pattern_forward = r"(?=([ACG][ACGT]{19}[ACGT]GG))"
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse = r"(?=(CC[ACGT][ACGT]{19}[TGC]))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"


# Defining the patterns used for secondary structures
pattern_RNAstructure = r".{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)"
pattern_RNAenergy = r"\s\((.+)\)"


# Thresholds used when processing secondary structures
low_energy_threshold = -30
high_energy_threshold = -18


# Threshold when looking at the off-target sites
offtarget_threshold = 75


# guide RNA
guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"


work_dir = os.environ['WORK_DIR']
dir = os.path.join(work_dir, "data/%s/" % args.gene)
dir_seq = os.path.join(dir, "outputs", "mm10db/")
dir_trg = dir_seq
dir_list = dir_seq


inputFile_exon = "%s_exon_sequences.txt" % args.chr
inputFile_chr = "%s.txt" % args.chr
outputFile = "potentialTargets.txt"
out_RNAfold = "RNAfold_output.txt"
out_targetsToScore = "targetsToScore.txt"
in_targetScores = "targets_scored.txt"
C_program = "./findMismatches_threads"
bowtie_index = os.path.join(dir, "indexes", "bowtie2", "%s_%s" % (args.gene, args.chr))

accepted_targets = dir_trg+"%s_accepted_targets.txt" % args.chr
accepted_targets_sorted = dir_trg+"%s_accepted_targets_sortedByGene.txt" % args.chr
accepted_targets_Excel = dir_trg+"%s_accepted_targets_ExcelFriendly.tsv" % args.chr
rejected_targets = dir_trg+"%s_rejected_targets.txt" % args.chr

exon_list = "%s_exonList.txt" % args.chr
offTargetSites = dir_seq+"%s_offtargetSites.txt" % args.chr
nb_threads_C = args.tc
nb_threads_Bowtie = args.tb

# This line has been commented out since it is not relevant for artificial
# genomes.
# chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "M", "X", "Y"]


padding = 17 # digestion site must be on exon, but the whole target does not have to.
# Note: if this value is different than that used to create exon sequences, we might have an error with trying to locate targets.


# Temporary files used in the script
tempTargetFile = "reads.txt"
alignmentFile = "alignedReads.txt"


# Parsing the arguments
# Given the use of argparse with default arguments, this parsing process is no
# longer necessary, so it was commented out.
"""
if len(sys.argv) != 4:
    print "Wrong number of arguments: 3 expected, "+str(len(sys.argv))+" given."
    print "Usage: python target_identitification.py nb_threads_C=<int> nb_threads_Bowtie=<int> genes=<all or gene name or set name>"
    quit ()
else:
    if "nb_threads_C=" in sys.argv[1]:
        nb_threads_C = sys.argv[1][13:]
    else:
        print "Error in first argument"
        print "Usage: python target_identitification.py nb_threads_C=<int> nb_threads_Bowtie=<int> genes=<all or gene name or set name>"
        quit ()
    if "nb_threads_Bowtie=" in sys.argv[2]:
        nb_threads_Bowtie = sys.argv[2][18:]
    else:
        print "Error in first argument"
        print "Usage: python target_identitification.py nb_threads_C=<int> nb_threads_Bowtie=<int> genes=<all or gene name or set name>"
        quit ()
    if "genes=" in sys.argv[3]:
            gene = sys.argv[3][6:]
            # If a gene (or gene set) is given has an argument, we only run the detection for that gene
"""
if args.list != "all":
    # renaming the input files
    inputFile_exon = inputFile_exon[:-4]+"_"+gene+".txt"
    exon_list = exon_list[:-4]+"_"+gene+".txt"
    # renaming the ouput files
    outputFile = outputFile[:-4]+"_"+gene+".txt"
    accepted_targets = accepted_targets[:-4]+"_"+gene+".txt"
    accepted_targets_sorted = accepted_targets_sorted[:-4]+"_"+gene+".txt"
    accepted_targets_Excel = accepted_targets_Excel[:-4]+"_"+gene+".tsv"
    rejected_targets = rejected_targets[:-4]+"_"+gene+".txt"

    print "Working with "+inputFile_exon+" and "+exon_list
    if os.path.isfile(dir_seq+inputFile_exon)==False:
        print "Error: file "+inputFile_exon+" does not exist."
        quit()
    if os.path.isfile(dir_list+exon_list)==False:
        print "Error: file "+exon_list+" does not exist."
        quit()
"""            
    else:
        print "Running the method on the whole genome."
        print "WARNING: the program might crash if the number of potential targets exceeds the memory available on this computer."
else:
    print "Error in third argument"
    print "Usage: python target_identitification.py nb_threads_C=<int> nb_threads_Bowtie=<int> genes=<all or gene name or set name>"
    quit ()
"""



#############################
##   Auxiliary functions   ##
#############################


# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq




# Function that replaces U with T in the sequence (to go back from RNA to DNA)
def transToDNA(rna):
    switch_UT = string.maketrans('U', 'T')
    dna = rna.translate(switch_UT)
    return dna




# Function that calculates the AT% of a given sequence
def AT_percentage(seq):
    total = 0.0
    length = float(len(seq))
    for c in seq:
        if c in "AT":
            total += 1
    return 100*total/length




# Function that aligns two given sequences and returns their similarity
def NeedlemanWunsch(seq1,seq2):
    d = -5 # Gap penalty
    A = seq1 # First sequence to be compared
    B = seq2 # Second ""
    I = range(len(seq1)) # To help iterate (Pythonic)
    J = range(len(seq2)) # ""
    F = [[0 for i in seq1] for j in seq2] # Fill a 2D array with zeroes
    # Similarity matrix from Wikipedia:
    S = \
    {'A': {'A': 10, 'G': -1, 'C': -3, 'T': -4},
        'G': {'A': -1, 'G':  7, 'C': -5, 'T': -3},
        'C': {'A': -3, 'G': -5, 'C':  9, 'T':  0},
        'T': {'A': -4, 'G': -3, 'C':  0, 'T':  8}}

    # Initialization
    for i in I:
        F[i][0] = d * i
    for j in J:
        F[0][j] = d * j

    # Scoring
    for i in I[1:]:
        for j in J[1:]:
            Match = F[i-1][j-1] + S[A[i]][B[j]]
            Delete = F[i-1][j] + d
            Insert = F[i][j-1] + d
            F[i][j] = max(Match, Insert, Delete)

    # Traceback
    AlignmentA = ""
    AlignmentB = ""
    i = len(seq1) - 1
    j = len(seq2) - 1

    while (i > 0 and j > 0):
        Score = F[i][j]
        ScoreDiag = F[i - 1][j - 1]
        ScoreUp = F[i][j - 1]
        ScoreLeft = F[i - 1][j]
        if (Score == ScoreDiag + S[A[i]][B[j]]):
            AlignmentA = A[i] + AlignmentA
            AlignmentB = B[j] + AlignmentB
            i -= 1
            j -= 1
        elif (Score == ScoreLeft + d):
            AlignmentA = A[i] + AlignmentA
            AlignmentB = "-" + AlignmentB
            i -= 1
        elif (Score == ScoreUp + d):
            AlignmentA = "-" + AlignmentA
            AlignmentB = B[j] + AlignmentB
            j -= 1
        else:
            print("algorithm error?")
    while (i > 0):
        AlignmentA = A[i] + AlignmentA
        AlignmentB = "-" + AlignmentB
        i -= 1
    while (j > 0):
        AlignmentA = "-" + AlignmentA
        AlignmentB = B[j] + AlignmentB
        j -= 1

    # Similarity
    lenA = len(AlignmentA)
    lenB = len(AlignmentB)
    sim1 = ""
    sim2 = ""
    len0 = 0
    k = 0
    total = 0.0
    similarity = 0.0

    if (lenA > lenB):
        sim1 = AlignmentA
        sim2 = AlignmentB
        len0 = lenA
    else:
        sim1 = AlignmentB
        sim2 = AlignmentA
        len0 = lenB

    while (k < len0):
        if (sim1[k] == sim2[k]):
            total += 1
        k += 1

    similarity = total / len0 * 100

#    print AlignmentA
#    print AlignmentB
#    print "\t\t\t"+str(similarity)
    return similarity




#######################
##   Main function   ##
#######################




###################################
##   Processing the input file   ##
###################################

print strftime("%H:%M:%S", localtime())+":\tGetting ready to process "+inputFile_exon
inFile = open(dir_seq+inputFile_exon,'r')

exonLineNumber=0
possibleTargets=dict()
removedTargets=dict()

# For every line in the input file

for line_exon in inFile:
    exonLineNumber+=1
    #print exonLineNumber
    
    # we parse the line and look for forward sequences
    match_exon = re.findall(pattern_forward,line_exon)
    if match_exon:
        for i in range(0,len(match_exon)):
            target23 = match_exon[i]
            if target23 in possibleTargets:
                possibleTargets[target23].append(exonLineNumber)
            else:
                possibleTargets[target23]=[]
                possibleTargets[target23].append(exonLineNumber)

    # we parse the line and look for reverse sequences
    match_exon = re.findall(pattern_reverse,line_exon)
    if match_exon:
        for i in range(0,len(match_exon)):
            target23 = rc(match_exon[i])
            if target23 in possibleTargets:
                possibleTargets[target23].append(exonLineNumber)
            else:
                possibleTargets[target23]=[]
                possibleTargets[target23].append(exonLineNumber)

inFile.close()
#print "\t\tSKIPPED THE IDENTIFICATION OF REVERSE SEQUENCES!"

print "\n"+strftime("%H:%M:%S", localtime())+":\t%d potential targets have been identified." % (len(possibleTargets))



##############################################################
##   Removing targets that have multiple matches in exons   ##
##############################################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that have been observed more than once."

targetsToRemove=[]
for target23 in possibleTargets:
    # number of occurrences of the target
    total_occurrences = len(possibleTargets[target23])
    
    # number of occurrences of the reverse complement target
    reverse_target23 = rc(target23)
    reverse_also_exists = False
    if reverse_target23 in possibleTargets:
        total_occurrences += len(possibleTargets[reverse_target23])
        reverse_also_exists = True
    
    # we reject if the total is greater than 1
    if total_occurrences>1:
        targetsToRemove.append(target23)
        # we also reject the reverse complement if it exists
        if reverse_also_exists:
            targetsToRemove.append(reverse_target23)

for target23 in targetsToRemove:
    # if the target is not already removed (as reverse-complement of another one)...
    if target23 in possibleTargets:
        # ... then we remove it
        del possibleTargets[target23]
        removedTargets[target23] = "Multiple matches in exons"

print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))



###############################################
##   Using Bowtie to find multiple matches   ##
###############################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tPreparing file for Bowtie analysis."
outFile = open(tempTargetFile,'w')

tempTargetDict_offset = dict()
for target23 in possibleTargets:
    similarTargets = [target23[0:20]+"AGG", target23[0:20]+"CGG", target23[0:20]+"GGG", target23[0:20]+"TGG", target23[0:20]+"AAG", target23[0:20]+"CAG", target23[0:20]+"GAG", target23[0:20]+"TAG"]
    for seq in similarTargets:
        outFile.write(seq+"\n")
        tempTargetDict_offset[seq] = target23
outFile.close()



print "\n"+strftime("%H:%M:%S", localtime())+":\tFile ready. Calling Bowtie."

cmd = "bowtie2 -x %s -p %d --reorder --no-hd -t -r -U %s -S %s" % (bowtie_index, nb_threads_Bowtie, tempTargetFile, alignmentFile)
call([cmd],shell=True)



print "\n"+strftime("%H:%M:%S", localtime())+":\tStarting to process the Bowtie results."
inFile = open(alignmentFile,'r')
bowtieLines = inFile.readlines()
inFile.close()

targetsToRemove=[]

i=0
while i<len(bowtieLines):
    nb_occurences = 0
    # we extract the read and use the dictionnary to find the corresponding target
    read = bowtieLines[i].rstrip().split("\t")[9]
    seq = ""
    if read in tempTargetDict_offset:
        seq = tempTargetDict_offset[read]
    elif rc(read) in tempTargetDict_offset:
        seq = tempTargetDict_offset[rc(read)]
    else:
        print "Problem? "+read


    # we count how many of the eight reads for this target have a perfect alignment
    for j in range(i,i+8):
        if "XM:i:0" in bowtieLines[j]:
            nb_occurences += 1
            # we also check whether this perfect alignment also happens elsewhere
            if "XS:i:0"  in bowtieLines[j]:
                nb_occurences += 1

    # if that number is at least two, the target is removed
    if nb_occurences > 1:
        targetsToRemove.append(seq)

    # we continue with the next target
    i+=8


# we can remove the dictionnary
del tempTargetDict_offset

for target23 in targetsToRemove:
    # if the target is not already removed (as reverse-complement of another one)...
    if target23 in possibleTargets:
        # ... then we remove it
        del possibleTargets[target23]
        removedTargets[target23] = "Multiple matches in genome"

print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))



##############################################
##   Removing targets that have AT% < 45%   ##
##############################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that have AT% strictly below 45%."

targetsToRemove=[]
for target23 in possibleTargets:
    target = target23[0:20]
    if AT_percentage(target)<45:
        targetsToRemove.append(target23)

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "AT%"


print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))



############################################
##   Removing targets that contain TTTT   ##
############################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that contain TTTT."

targetsToRemove=[]
for target23 in possibleTargets:
    if "TTTT" in target23:
        targetsToRemove.append(target23)

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "TTTT"


print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))



############################################
##   Removing targets that contain TTT   ##
############################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that contain TTT."

print "\t\tSKIPPED!"
#targetsToRemove=[]
#for target23 in possibleTargets:
#    if "TTT" in target23:
#        targetsToRemove.append(target23)
#
#for target23 in targetsToRemove:
#    del possibleTargets[target23]
#    removedTargets[target23] = "TTT"


print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))



################################################################
##   Removing targets that are too close to reverser primer   ##
################################################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that are too close to reverse primer."

#print "\t\tSKIPPED!"
targetsToRemove=[]
for target23 in possibleTargets:
    target = target23[0:20]
    if NeedlemanWunsch(target,"AAAAGCACCGACTCGGTGCC")>60:
        targetsToRemove.append(target23)

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "Too close to reverse primer"


print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))



##########################################
##   Calculating secondary structures   ##
##########################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tCalculating secondary structures."

# WARNING: removing any existing version of the RNAfold output file.
call(["rm -f "+out_RNAfold],shell=True)

# Calling RNAfold on all targets

temp_counter = 0
for target23 in possibleTargets:
    target = "G"+target23[1:20]
    structure = target+guide
    cmd = "echo "+structure+" | RNAfold --noPS >> "+out_RNAfold
    call([cmd],shell=True)
    temp_counter+=1
    if (temp_counter%100000)==0:
        print strftime("%H:%M:%S", localtime())+":\t\t"+str(temp_counter)+" targets processed."
#    if temp_counter>break_point:
#        print strftime("%H:%M:%S", localtime())+":\t\tReaching a break point!"
#        break

total_number_structures = temp_counter


#########################################
##   Processing secondary structures   ##
#########################################

print "\n"+strftime("%H:%M:%S", localtime())+":\tProcessing secondary structures."

inFile = open(out_RNAfold,'r')
RNA_structures = inFile.readlines()
inFile.close()

targetsToRemove=[]
i=0
for target23 in possibleTargets:
    L1 = RNA_structures[2*i].rstrip()
    L2 = RNA_structures[2*i+1].rstrip()
    target = L1[:20]
    if transToDNA(target) != target23[0:20] and transToDNA("C"+target[1:]) != target23[0:20] and transToDNA("A"+target[1:]) != target23[0:20]:
        print "Error? "+target23+"\t"+target
        quit()
#    print L1
#    print target
#    print L2
    match_structure = re.search(pattern_RNAstructure,L2)
    if match_structure:
        # The structure is correct, we only reject if the energy is too low
        energy = ast.literal_eval(match_structure.group(1))
        if energy < low_energy_threshold:
            targetsToRemove.append(transToDNA(target23))
    else:
        match_energy = re.search(pattern_RNAenergy,L2)
        if match_energy:
            # The structure is not correct, we only reject if the energy is not high enough
            energy = ast.literal_eval(match_energy.group(1))
            if energy <= high_energy_threshold:
                targetsToRemove.append(transToDNA(target23))
    i+=1
#    if i>break_point-1:
#        print strftime("%H:%M:%S", localtime())+":\t\tReaching a break point!"
#        break

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "Secondary structure or energy"


print "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))




#########################
##   Scoring targets   ##
#########################

print "\n"+strftime("%H:%M:%S", localtime())+":\tScoring targets."

i=0
targetsToRemove=[]

print "\n"+strftime("%H:%M:%S", localtime())+":\tWriting targets to file."

outFile = open(out_targetsToScore,'w')

temp_counter = 0
for target23 in possibleTargets:
#    print strftime("%H:%M:%S", localtime())+":\t\t"+target23
    target = target23[0:20]
    outFile.write(target+"\n")
    temp_counter += 1
outFile.close()

total_number_scores = temp_counter

print "\n"+strftime("%H:%M:%S", localtime())+":\tCalling C program.\n"

cmd = "%s %d %s %s %s %d" % (C_program, nb_threads_C, out_targetsToScore, offTargetSites, in_targetScores, offtarget_threshold)
stderr_fd = open(os.path.join(dir_seq, "stderr.txt"), 'w')
call([cmd],shell=True, stderr=stderr_fd)
stderr_fd.close()

print "\n"+strftime("%H:%M:%S", localtime())+":\tReading the results from file, and processing."

inFile = open(in_targetScores,'r')
outFile = open(tempTargetFile,'w')

for target23 in possibleTargets:
    line = inFile.readline().rstrip()
    t20 = line.split("\t")[0]
    score = ast.literal_eval(line.split("\t")[1])
    if t20 != target23[0:20]:
        print "Problem? "+t20+" - "+target
        quit()
#    print target23+"\t"+str(score)
    if score <offtarget_threshold:
        targetsToRemove.append(target23)
    else:
        outFile.write(target23+"\n")

inFile.close()
outFile.close()


for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "Off-target score"


print "\t\t%d potential targets are selected as successful candidates." % (len(possibleTargets))


print "\n"+strftime("%H:%M:%S", localtime())+":\tCalculating their exact position using Bowtie."

cmd = "bowtie2 -x %s -p %d --reorder --no-hd -t -r -U %s -S %s" % (bowtie_index, nb_threads_Bowtie, tempTargetFile, alignmentFile)
call([cmd],shell=True)



print "\n"+strftime("%H:%M:%S", localtime())+":\tSaving the results."

outFile = open(accepted_targets,'w')
inFile_score = open(in_targetScores,'r')
inFile_Bowtie = open(alignmentFile,'r')

inFile = open(dir_list+exon_list,'r')
exons = inFile.readlines()
inFile.close()


i=0
for target23 in possibleTargets:

    # For each target, we save...

    output_line = ""
    target = target23[0:20]
    
    # ... the sequence
    output_line += target+"\t"
    
    # ... the secondary and the energy
    while i<total_number_structures:
        L1 = RNA_structures[2*i].rstrip()
        L2 = RNA_structures[2*i+1].rstrip()
        target_RNA = L1[:20]
        if transToDNA(target_RNA) == target or transToDNA("C"+target_RNA[1:]) == target or transToDNA("A"+target_RNA[1:]) == target:
            structure = L2.split(" ")[0]
            energy = L2.split(" ")[1][1:-1]
            output_line += L1+"\t"+structure+"\t"+energy+"\t"
            break
        i+=1
    if i == total_number_structures:
        print "Error? "+target+" not found in "+out_RNAfold
        quit()

    # ... the off-target score
    j=0
    while j<total_number_scores:
        line = inFile_score.readline().rstrip()
        t20 = line.split("\t")[0]
        score = line.split("\t")[1]
        if t20 == target:
            output_line += score+"\t"
            break
        j+=1
    if j == total_number_scores:
        print "Error? "+target+" not found in "+in_targetScores
        quit()

    # ... the position
    line = inFile_Bowtie.readline().rstrip()
    tempArray = line.split("\t")
    chr = tempArray[2]
    pos = ast.literal_eval(tempArray[3])
    seq = tempArray[9]
    match = tempArray[10]
    if match != "IIIIIIIIIIIIIIIIIIIIIII":
        print "Error? Imperfect Bowtie match for target "+target23
        quit()
    if seq == target23:
        output_line += chr+"\t"+str(pos)+"\t"+str(pos+19)+"\t+\t"
    elif rc(seq) == target23:
        output_line += chr+"\t"+str(pos+3)+"\t"+str(pos+22)+"\t-\t"
    else:
        print "Error? "+target23+" not found in Bowtie results."
        print "Seq:\t"+seq
        print "Target23:\t"+target23
        quit()

    # ... the gene name, and position respective to the CDS
    j=0
    while j<len(exons):
        tempArray = exons[j].rstrip().split("\t")
        chr_exon = tempArray[1]
        start_exon = ast.literal_eval(tempArray[3])
        end_exon = ast.literal_eval(tempArray[4])
        if chr == chr_exon and start_exon-padding <= pos and pos <= end_exon+padding:
            # gene name
            output_line += tempArray[0]+"\t"
            # position w.r.t. CDS
            CDS_start = ast.literal_eval(tempArray[5])
            CDS_end = ast.literal_eval(tempArray[6])
            CDS_flag = tempArray[7]
            gene_direction = tempArray[2]
            if CDS_start  <= pos and pos <= CDS_end:
                output_line += "CDS\t"
            else:
                if CDS_flag == "CDS_Intersection":
                    output_line += "Isoform-dependent\t"
                else:
                    if pos < CDS_start:
                        if gene_direction=="+":
                            output_line += "5'-UTR\t"
                        else:
                            output_line += "3'-UTR\t"
                    else:
                        if gene_direction=="+":
                            output_line += "3'-UTR\t"
                        else:
                            output_line += "5'-UTR\t"
            break
        j+=1
    if j == len(exons):
        print "Error? "+target+" not found in "+exon_list
        quit()

    # ... and the corresponding primers
    forward_primer_1 = "CACTATAGG"+target[1:]+"gttttagagctaGAAAtagc"
    forward_primer_2 = "gggccTAATACGACTCACTATAGG"+target[1:]+"g"
    output_line += forward_primer_1+"\t"+forward_primer_2+"\n"

    outFile.write(output_line)

inFile_score.close()
inFile_Bowtie.close()
outFile.close()

outFile = open(rejected_targets,'w')
for target23 in removedTargets:
    outFile.write(target23+"\t"+removedTargets[target23]+"\n")
outFile.close()

print "\n"+strftime("%H:%M:%S", localtime())+":\tPreparing the Excel-friendly output."

cmd = "sort -k10,10 -k7,7n "+accepted_targets+" > "+accepted_targets_sorted
call([cmd],shell=True)

inFile = open(accepted_targets_sorted,'r')
outFile = open(accepted_targets_Excel,'w')

for line in inFile:
    tempArray = line.rstrip().split("\t")
    outFile.write(tempArray[0]+"\t"+tempArray[1]+"\t"+tempArray[3]+"\t"+tempArray[4]+"\t"+tempArray[5]+":"+tempArray[6]+"-"+tempArray[7]+"\t"+tempArray[8]+"\t"+tempArray[9]+"\t"+tempArray[10]+"\t"+tempArray[11]+"\t"+tempArray[12]+"\n\t"+tempArray[2]+"\n")

inFile.close()
outFile.close()

print "\n"+strftime("%H:%M:%S", localtime())+":\tDone."


#######################
##                   ##
##    End of File    ##
##                   ##
#######################
