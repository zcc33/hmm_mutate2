#!/usr/bin/env python
#Mutates a protein sequence according to state and transition probabilities of a representative HMM
#Input1: .hmm file
#Input2: sequence file (in FASTA format)
#Input3: percent similarity - should be, .e.g., 90 for 90% similarity (which implies 10% mutation)
#Input4: number of mutated sequences to generate
#Algorithm shaves off overhanging deletions at the beginning/end;
#Mutates upon a single match character, a block of consecutive insertions, or a block of consecutive deletions
#By referring to the HMM's previous Match state, transitioning and generating one character (be it match, insertion, or deletion)
###################################################
import os
import sys
import math
import random

#Fill in location of hmmsearch program in HMMER 3.1b1 or 3.1b2
#hmm_search = "/Users/yancyliao/Dropbox/bioinformatics/hmmer-3.1b2-macosx-intel/src/hmmsearch"
hmm_search = "/fs/sw/RedHat-7-x86_64/common/local/hmmer/3.1b1/bin/hmmsearch"

hmm_file = sys.argv[1]
seq_file = sys.argv[2]
pct_sim = float(sys.argv[3])
num_seq = int(sys.argv[4])
output_file = "hmmsearch_out"

pct_change = (100-pct_sim)*0.01
new_seqs = ["" for i in range(num_seq)]
    
#some variable initializations
amino = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
length = 0
sequence = ""
identifier = ""

#strip "-" from sequence
with open(seq_file, "r") as sequence_file:
    old = sequence_file.read()

old_text = old.split("\n")[0] + "\n"
for i in range(1, len(old.split("\n"))):
    old_text = old_text + old.split("\n")[i].replace("-", "") + "\n"

seq_file = seq_file + "~"
with open(seq_file, "w") as sequence_file:
    sequence_file.write(old_text)
    
#run hmm search
os.system(hmm_search + " " + hmm_file + " " + seq_file + " > " + output_file)

#get the sequence and identifier
with open(seq_file, "r") as sequence_file:
    for line in sequence_file:
        if not line[0] == ">":
            sequence = sequence + line
        else:
            identifier = line
sequence = sequence.replace("\n", "")
identifier = identifier.replace("\n", "").replace(">", "")

#split the identifier (must be in the form A|B, A|B|C, or A|B/C)
split = identifier.split("|")
if len(split) == 1:
    raise ValueError('Identifier must be in the form A|B, A|B|C, or A|B/C')
elif len(split) ==2:
    split2 = identifier.split("/")
    if len(split2) == 1:
        id1 = identifier.split("|")[0] + "|"
        id2 = identifier.split("|")[1]
        id3 = ""
    else:
        id1 = identifier.split("|")[0] + "|"
        id2 = identifier.split("/")[0].split("|")[1]
        id3 = "/" + identifier.split("/")[1]
else:
    id1 =  identifier.split("|")[0] + "|"
    id2 = identifier.split("|")[1]
    id3 = ""
    for i in range(2, len(split)):
        id3 = id3 + "|" + split[i]

#name the output file
output_file2 = "mutated_" + id2 + "_" + str(int(pct_sim))+"similarity.fa"

#get the length of the hmm
with open(hmm_file, "r") as data:
    for line in data:
        if line.split()[0] == "LENG":
            length = int(line.split()[1])


#collect all the probability information from the hmm
#note: matrix[position][amino acid]
mmatrix = [[0 for x in range(20)] for y in range(length)]
imatrix = [[0 for x in range(20)] for y in range(length)]
transitions = [[0 for x in range(7)] for y in range(length)]
start_insertions = [0 for x in range(20)]
start_transitions = [0 for x in range(7)]

with open(hmm_file, "r") as data:
    index = 1
    lines = data.readlines()
    
    for j in range(len(lines)):
        
        if lines[j].split()[0] == "COMPO":
            for i in range(20):
                start_insertions[i] = math.exp(-float(lines[j+1].split()[i]))
            for i in range(7):
                if lines[j+2].split()[i] == "*":
                    start_transitions[i] = "*"
                else:
                    start_transitions[i] = math.exp(-float(lines[j+2].split()[i]))
        
        if lines[j].split()[0] == str(index):
            
            #fill in the match matrix
            for i in range(20):
                mmatrix[index-1][i] = math.exp(-float(lines[j].split()[i+1]))
                
            #fill in the insertion matrix
            for i in range(20):
                imatrix[index-1][i] = math.exp(-float(lines[j+1].split()[i]))
                
            #fill in the transition matrix
            for i in range(7):
                if lines[j+2].split()[i] == "*":
                    transitions[index-1][i] = "*"
                else:
                    transitions[index-1][i] = math.exp(-float(lines[j+2].split()[i]))
                    
            index = index+1
            
        if index > length:
            break
        
        
#parse the output file
first_hmm_pos = 0
first_seq_pos = 0
last_hmm_pos = 0
last_seq_pos = 0

with open(output_file, "r") as out:
    
    lines = out.readlines()
    first_ind = True
    traversal = ""
    hmm_traversal = ""
    
    for i in range(len(lines)):
        if len(lines[i]) > 1:
            if len(lines[i].split()[0]) > 5 and identifier.startswith(lines[i].split()[0]):
                traversal += lines[i].split()[2]
                hmm_traversal+= lines[i-2].split()[2]
                if first_ind:
                    first_hmm_pos = int(lines[i-2].split()[1])
                    first_seq_pos = int(lines[i].split()[1])
                    
                    first_ind = False
                last_hmm_pos = int(lines[i-2].split()[3])
                last_seq_pos = int(lines[i].split()[3])
 
    #use hmm_traversal to get rid of filler that was put into sequence and traversal that don't correspond to any hmm position
    #for i in range(len(hmm_traversal)):
    #    if hmm_traversal[i] == "." and not traversal[i] == "-":
    #        traversal = traversal[:i] + traversal[i+1:]
    #        sequence = sequence[:i+first_seq_pos-1] + sequence[i+first_seq_pos:]

#generate each of the required sequences via mutation
#the algorithm
for f in range(num_seq):

    #process the stuff above
    characters = ["" for i in range(length)]
    states = ["" for i in range(length)]
    hmm_index = [i for i in range(length)]
        
    #fill in the beginning truncated
    for i in range(1, first_hmm_pos):
        if abs(i - first_hmm_pos) > abs(first_seq_pos - 1):
            characters[i-1] = sequence[i-1]
            states[i-1] = "deletion"
        else:
            characters[i-1] = sequence[first_seq_pos - abs(i-first_hmm_pos)-1]
            states[i-1] = "match"
            
    #fill in the end truncated
    for i in range(last_hmm_pos+1, length+1):
        if(i - last_hmm_pos > len(sequence) - last_seq_pos):
            characters[i-1]= "-"
            states[i-1] = "deletion"
        else:
            characters[i-1] = sequence[last_seq_pos + (i - last_hmm_pos) - 1]
            states[i-1] = "match"

                    
    #fill in the middle traversed parts, collapsing consectuvive deletions and insertions w/r/t hmm
    t_index = 0
    i = first_hmm_pos
    while t_index < len(traversal):
        if traversal[t_index] == "-":
            characters[i-1] += "-"
            states[i-1] = "deletion"
            while traversal[t_index+1] == "-":
                characters = characters[:i] + characters[i+1:]
                states = states[:i] + states[i+1:]
                hmm_index = hmm_index[:i] + hmm_index[i+1:]
                characters[i-1]+="-"
                t_index+=1
        elif traversal[t_index].isupper():
            characters[i-1] += traversal[t_index]
            states[i-1] = "match"
        else:
            characters = characters[:i-1] + [traversal[t_index]] + characters[i-1:]
            states = states[:i-1] + ["insertion"] + states[i-1:]
            hmm_index = hmm_index[:i-1] + [i-2] + hmm_index[i-1:]
            while traversal[t_index+1].islower():
                characters[i-1] += traversal[t_index+1]
                t_index+=1
        i+=1
        t_index+=1
    
    #randomly pick positions of the character array to perform mutations on, simulating a match state in the previous position
    req_changes = int(length*pct_change)
    changed = [False for i in range(len(characters))]
    while(sum(changed) < req_changes):
        pick = random.randint(0, len(characters)-1)
        index = hmm_index[pick]
        if not changed[pick]:
            if pick ==0:
                #let's take a look at the start position
                roll = random.uniform(0,1)
                dummy = 0
                choice = -1
                #decide whether to go to a match, insertion, or deletion
                for i in range(3):
                    dummy += start_transitions[i]
                    if roll < dummy:
                        choice = i
                        break
                #match
                if choice == 0:
                    roll = random.uniform(0,1)
                    dummy = 0
                    match_char = ""
                    for j in range(20):
                        dummy+= mmatrix[index][j]
                        if roll < dummy:
                            match_char = amino[j]
                            break
                    if not (characters[pick] == match_char and states[pick] == "match"):
                        changed[pick] = True
                    states[pick] = "match"
                    characters[pick] = match_char
                    
                #insertion
                elif choice == 1:
                    roll = random.uniform(0,1)
                    dummy = 0
                    insertion_char = ""
                    for j in range(20):
                        dummy+= start_insertions[j]
                        if roll < dummy:
                            insertion_char = amino[j]
                            break
                    if not (characters[pick] == insertion_char and states[pick] == "insertion"):
                        changed[pick] = True
                    states[pick] = "insertion"
                    characters[pick] = characters[pick] + insertion_char
                    
                elif choice == 2:
                    if not states[pick] == "deletion":
                        changed[pick] = True
                    states[pick] = "deletion"
                    characters[pick] == "-"
                else:
                    sys.exit("Random roll did not work")
                    
                    
            else:
                #let's take a look at the pick-1 position
                roll = random.uniform(0,1)
                dummy = 0
                choice = -1
                #decide whether to go to a match, insertion, or deletion
                for i in range(3):
                    dummy += transitions[index-1][i]
                    if roll < dummy:
                        choice = i
                        break
                #go to match
                if choice ==0:
                    roll = random.uniform(0,1)
                    match_char = ""
                    dummy = 0
                    for j in range(20):
                        dummy+= mmatrix[index][j]
                        if roll < dummy:
                            match_char = amino[j]
                            break
                    if not (characters[pick] == match_char and states[pick] == "match"):
                        changed[pick] = True
                    states[pick] = "match"
                    characters[pick] = match_char
                #go to insert
                elif choice==1:
                    roll = random.uniform(0,1)
                    dummy = 0
                    insertion_char = ""
                    for j in range(20):
                        dummy+= imatrix[index][j]
                        if roll < dummy:
                            insertion_char = amino[j]
                            break
                    if not (characters[pick] == insertion_char and states[pick] == "insertion"):
                        changed[pick] = True
                    states[pick] = "insertion"
                    characters[pick] = characters[pick] + insertion_char
                #go to delete
                else:
                    if not states[pick] == "deletion":
                        changed[pick] = True
                    states[pick] = "deletion"
                    characters[pick] == "-"
                            
    new_seq = ""
    for i in range(len(characters)):
        new_seq += characters[i]
    
    new_seqs[f] = new_seq.replace("-", "").upper()

#output the newly generated sequences in .fa output file
file = open(output_file2, "w")
for i in range(num_seq):
    file.write(">" + id1 + id2 + "_m" + str(i+1) +id3 +"\n")
    len_remain = len(new_seqs[i])
    index = 0
    while(len_remain > 80):
        file.write(new_seqs[i][index:index+80]+"\n")
        index+=80
        len_remain-=80
    file.write(new_seqs[i][index:]+"\n")
file.close()


print "HMM Mutate has run successfully.\n" + str(num_seq) + " sequences of " + str(sys.argv[3]) + "% similarity have been saved in " + output_file2 +"."