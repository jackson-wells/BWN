import sys
from Bio.SubsMat import MatrixInfo

def printMatrixFile(matrixName):
    OUTPUT = open(matrixName+".txt",'w')
    S = getMatrix(matrixName)
    AA =['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
#    print >>OUTPUT, "#",
#    for a in AA:
#	print >>OUTPUT, "\t", a,
#    print >>OUTPUT, ""
    # rows of matrix
    for a in AA:
#        print >>OUTPUT, a,
        # columns of matrix
        for b in AA:
            print >>OUTPUT, "\t", getScore(S,a,b),
        print >>OUTPUT, ""
    
def getMatrix(matrixName):
    if matrixName == 'blosum62':
         S = MatrixInfo.blosum62
         return S
    elif matrixName == 'blosum90':
         S = MatrixInfo.blosum90
         return S
    elif matrixName == 'pam30':
         S = MatrixInfo.pam30
         return S
    elif matrixName == 'pam60':
         S = MatrixInfo.pam60
         return S
    elif matrixName == 'pam250':
         S = MatrixInfo.pam250
         return S
    else:
        print "Unknown matrix name. Supported names are:"
        print "blosum62"
        # add new matrix names as they become available


def getScore(S,a,b):
    if (a,b) not in S:
        return S[b,a]
    return S[a,b]





########
# MAIN #
########

usage = "usage: " + sys.argv[0] + " <matrix name>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

matrixName = sys.argv[1]

printMatrixFile(matrixName)
