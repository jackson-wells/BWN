import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from Bio.Seq import Seq
from Bio import SeqIO
import random
import numpy as np
import time

usage = "usage: " + sys.argv[0] + " <database FASTA> <number of substrings> <start length> <end length> <algorithm>"
if len(sys.argv) != 6:
	print usage
    	sys.exit()

inputFasta = sys.argv[1]
count = int(sys.argv[2])
startLength = int(sys.argv[3])
endLength = int(sys.argv[4])
algorithm = sys.argv[5]

sequencesList = []
sequences = SeqIO.parse(inputFasta,'fasta')
for record in sequences:
        sequencesList.append((record.id,record.seq))

pid = str(os.getpid())
fh = "midScoredAnalysisRuntime." + pid + ".txt"
f = open(fh,"a")

data = {} # dictionary of lists
for i in range(count):
        rID,rSeq = random.choice(sequencesList)
	seqMid = len(rSeq)/2
        for k in range(startLength,endLength+1):
                # substring of length k
                substring = rSeq[seqMid:seqMid+k]
                # run job and get runtime
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

plt.xlabel("Length of randomly selected query from index (aa)")
plt.ylabel("Program runtime using " + algorithm + " algorithm (s)")
figName = "midScored." + pid + ".pdf"
plt.savefig(figName)#, bbox_inches='tight') 

