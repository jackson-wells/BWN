import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import random
import numpy as np
import time

#usage = "usage: " + sys.argv[0] + " <database FASTA> <number of substrings> <start length> <end length>"
usage = "usage: " + sys.argv[0] + " <database FASTA> <number of substrings> <max lenght of gap seq> <seq lenght>"
if len(sys.argv) != 5:
	print usage
    	sys.exit()

inputFasta = sys.argv[1]
count = int(sys.argv[2])
numMM = int(sys.argv[3]) + 1
seqLength = int(sys.argv[4])

sequencesList = []
sequences = SeqIO.parse(inputFasta,'fasta')
for record in sequences:
        sequencesList.append((record.id,record.seq))



#gets list of amino acids
alphabet = IUPAC.IUPACProtein.letters.upper()	
#creates random AA seq of designated length
#gapString = ''.join(random.choice(alphabet) for i in range (10))




pid = str(os.getpid())
fh = "mmScoredAnalysisRuntime." + pid + ".txt"
f = open(fh,"a")

data = {} # dictionary of lists
for k in range(0,numMM):
	for i in range(count):

		#new randomly selected seq 0:count times until gap length is increased

	        rID,rSeq = random.choice(sequencesList)
	        rStart = random.randint(0,len(rSeq)-seqLength)
		substring = rSeq[rStart:seqLength+rStart]

                # run job and get runtime
		
		randGapSeq = ''.join(random.choice(alphabet) for j in range (k))	
	
		insertionLocation = random.randint(0,seqLength)
		substring = substring[:insertionLocation] + randGapSeq + substring[insertionLocation:]

                command = "bwp-search -S -s " + str(substring)
                startTime = time.time()
                os.system(command)
                runtime = time.time() - startTime
                if k not in data:
                        data[k] = []
                data[k].append(float(runtime))
		f.write(str(k) + " " + str(substring) + " " + str(runtime) + "\n")
                print k, substring, runtime
                
f.close()
print data
xVals = []
mean = []
stdErr = []
# initialize curves list
curves = {}
# collect data
for k in data:
        xVals.append(k)
        mean.append(np.mean(data[k]))
        stdErr.append(np.std(data[k]))
        for i in range(len(data[k])):
                if i not in curves:
                        curves[i] = []
                curves[i].append(data[k][i])
# plot mean, errorbars                
plt.errorbar(xVals, mean, yerr=stdErr, fmt='o', ms=3, elinewidth=1, color='steelblue')
# plot individual curves
for i in range(count):
        plt.plot(xVals,curves[i], linewidth=1, color='lightblue')

plt.xlabel("Length of single random gap within random query of length " + str(seqLength)  + " (aa)")
plt.ylabel("Program runtime using gap-aware algorithm (s)")
figName = "mmScored." + pid + ".pdf"
plt.savefig(figName)#, bbox_inches='tight') 

