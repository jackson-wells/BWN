import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from Bio.Seq import Seq
from Bio import SeqIO
import random
import numpy as np
import time

usage = "usage: " + sys.argv[0] + " <output file>"
if len(sys.argv) != 2:
	print usage
    	sys.exit()

inputFile = sys.argv[1]

inputFilePID = inputFile.split(".")
inputFilePID = inputFilePID[1]

data = {}

fh = open(inputFile,'r')

line = fh.readline()

lineList = line.split()

length = len(lineList[1])

sampleCount = 0

while line:
	k = int(lineList[0])
	if k not in data:
		data[k] = []
	data[k].append(float(lineList[2]))
	line = fh.readline()
	lineList = line.split()

fh.close()


                
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

for i in curves:
        plt.plot(xVals,curves[i], linewidth=1, color='lightblue')

plt.xlabel("Number of single character random gap within random query of length " + str(length)  + " (aa)")
plt.ylabel("Program runtime using edit-distance algorithm (s)")
figName = "edGapFF." + inputFilePID + ".pdf"
plt.savefig(figName)#, bbox_inches='tight') 

