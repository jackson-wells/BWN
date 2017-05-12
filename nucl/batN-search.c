#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "batN-search.h"

int MAX_LINE_LENGTH = 10000000;
#define INTERVAL_FILE "index.batN"
char OUTPUT_FILE[100];

/* MEMORY ALLOCATION FUNCTIONS */

struct query initializeInputStruct(int sCount, int *sLength)
{
        struct query input;
        int i;
        int n = 100;
        input.length = (int *) malloc(sCount * sizeof(int));
        input.name = (char **) malloc(sCount * sizeof(char *));
        input.sequence = (char **) malloc(sCount * sizeof(char *));
        for(i = 0; i < sCount;i++)
        {
                input.name[i] = (char *) malloc(n*(sizeof(char)));
                input.sequence[i] = (char *) malloc(sLength[i] * sizeof(char));
        }
        return input;
}

/* HELPER FUNCTIONS */

/*  seqCount
 *
 *      Purpose: Find the number of sequences in the input file
 *              Returns: integer
 *              */
int seqCount(char *fileName)    /*reutrns number of sequences present in a multi-fasta file*/
{
        FILE *file = fopen(fileName,"r");
        int lineCount = 0;
        char *temp = (char *) malloc((MAX_LINE_LENGTH+1) * sizeof(char));
        while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)
        {
                lineCount++;
        }
        fclose(file);
        int seqCount = lineCount/2;
        /*free(temp);*/
        return seqCount;
}

/* removePrefix
 *
 * 	Prupose: remove the first two characters from a string
 * 		Returns: input string without its first two chars
 * 		*/
char *removePrefix(char *input)
{
	memmove(input, input+2, strlen(input));
	return input;
}

/*  seqLength
 *
 *      Purpose: get the lengths of input sequences from file
 *              Returns: array of integers
 *              */
int *seqLength(char *fileName,int seqCount,int cCount) /*returns an array of sequence lengths*/
{
        FILE *file = fopen(fileName,"r");
        char *seq = (char *) malloc(cCount * sizeof(char));
        int *seqLength = (int *) malloc(seqCount * sizeof(int));;
        int i = 0;
        while(fgets(seq,cCount-1,file) != NULL)
        {
            if(seq[0] != '>' && seq[0] != '\n')
            {
                seqLength[i] = strlen(seq) + 1;
                i++;
            }
        }
        fclose(file);
        /*free(seq);*/
        return seqLength;
}

/*  read_fasta
 *
 *      Purpose: Store FASTA file information into memory
 *              Returns: Nothing, Variables passed by reference
 *              */
void read_fasta(char *fileName, struct query *input)
{
    int charCount = MAX_LINE_LENGTH;
    char *temp = (char *) malloc(charCount * sizeof(char));
    /*m = newArr(seqC,sLength);*/
    FILE *file = fopen(fileName,"r");
    int i = 0;
    while(fgets(temp,charCount,file) != NULL)   /*loops through each line of a file*/
    {
            if(temp[0] == '>')      /*if line is a header*/
            {
                    strtok(temp,"\n");
                    memmove(temp, temp+1, strlen(temp));
                    strcpy(input->name[i],temp);
            }
            else if(temp[0] == '\n') /*if line is empty*/
            {
                    continue;
            }
            else    /*if line contains a nucleotide sequence*/
            {
                    strtok(temp,"\n"); /*strings read from file have extra \n added by file read*/
                    strcat(temp,"$");
                    strcpy(input->sequence[i],temp);    /*saving string in memory*/
                    input->length[i] = strlen(temp);
                    i++;
            }
    }
    fclose(file);
        /*free(temp);*/
}

void handleS(struct query *input, char *sequence)
{
    printf("\nHere is the string you entered:\n%s\n\n",sequence);
    int *seqL = (int *) malloc(sizeof(int));
    seqL[0] = strlen(sequence) + 1;
    int seqC = 1;
    strcat(sequence,"$");
    *input = initializeInputStruct(seqC,seqL);
    input->length[0] = seqL[0];
    strcpy(input->sequence[0],sequence);
    strcpy(input->name[0],"Command Line String Input");
}

void handleF(struct query *input,char *fileName)
{
    	int sCount = seqCount(fileName);
    	int cCount = MAX_LINE_LENGTH;
    	int *sLength = seqLength(fileName,sCount,cCount);
    	*input = initializeInputStruct(sCount,sLength);
    	read_fasta(fileName, input);
}

struct query manageInputs(char *argv[], int argc, int *sCount)
{
	struct query input;
        int c;
        if(argc <= 1)   /*no arguements supplied*/
        {
                printf("Please provide the necessary options\n\nuse -h for usage statement\n");
                exit(0);
        }
        opterr = 0;
	while ((c = getopt (argc, argv, "hf:s:m:o:")) != -1) /*options must be added here to be recognized, options followed by : take in a parameter*/
        {
                switch (c)
                {
                        case 'h':
                                printf("\nBurrows Wheeler Protein Alligner\n\nUsage: \"batN-search <options>\"\n\nOptions:\n\n-f\t\tFor input of a fasta file\n-s\t\tFor input of a string\n-h\t\tFor this usage statement\n-m\t\tTo designate maximum sequence length according to character count\n-o\t\tSpecify output file name (exclude file extentions)\n\n");
                                exit(0);

                        case 'f':
                                handleF(&input,optarg);
                                *sCount = seqCount(optarg);
                                break;

                        case 's' :
                                handleS(&input,optarg);
                                *sCount = 1;
                                break;

                        case '?' :
                                if(optopt == 's')
                                {
                                    fprintf (stderr, "Option -%c requires an argument.\n\nuse -h for usage statement\n", optopt);
                                }
                                else if(optopt == 'f')
                                {
                                    fprintf (stderr, "Option -%c requires an argument.\n\nuse -h for usage statement\n", optopt);
                                }
                                else if(isprint (optopt))
                                {
                                    fprintf (stderr, "Unknown option `-%c'.\n\nuse -h for usage statement\n", optopt);
                                }
                                exit(0);
                        case 'm' :
                                MAX_LINE_LENGTH = atoi(optarg);
			case 'o' :
                                strcpy(OUTPUT_FILE,optarg);

                }
        }
        return input;
}

int getCount()
{
	FILE *file = fopen(INTERVAL_FILE,"r");
	char temp[10];
	fgets(temp,10,file);
	fclose(file);
	removePrefix(temp);
	return atoi(temp);
}

int *getLength(int seqCount)
{
	FILE *file = fopen(INTERVAL_FILE,"r");
        int *length = (int *) malloc(seqCount*sizeof(int));
	char *temp = (char *) malloc((MAX_LINE_LENGTH+1)*sizeof(char));
	int i = 0;
	while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)
	{
		if(temp[0] == 'l')
		{
			removePrefix(temp);
			length[i] = atoi(temp);
			i++;
		}
	}
        fclose(file);
        return length;
}

struct index *getIndex(int *IseqCount)
{
	int seqCount = getCount();
	int *seqLength = getLength(seqCount);
	*IseqCount = seqCount;
	int charCount = MAX_LINE_LENGTH;
	char *temp = (char *) malloc((MAX_LINE_LENGTH+1)*sizeof(char));
	struct index *tempIndex = (struct index *) malloc(seqCount * sizeof(struct index));
	FILE *file = fopen(INTERVAL_FILE,"r");
    	int i,j;
	for(i = 0; i < seqCount; i++)
	{
		tempIndex[i].O = (int **) malloc(5 * sizeof(int *));
		for(j = 0; j < 5; j++)
		{
			tempIndex[i].O[j] = (int *) malloc(seqLength[i] * sizeof(int));
		}
		tempIndex[i].length = seqLength[i];
	}
	i = 0;
	int oCount = 0;
    	while(fgets(temp,charCount,file) != NULL)   /*loops through each line of a file*/
    	{
		strtok(temp,"\n");
		if(temp[0] == '\n')
		{
			continue;
		}
		else if(temp[0] == 'n')
		{
			continue;
		}
		else if(temp[0] == 'd')
		{
			removePrefix(temp);
			tempIndex[i].desc = (char *) malloc(strlen(temp)*sizeof(char));
			strcpy(tempIndex[i].desc,temp);
		}
		else if(temp[0] == 'l')
                {
			continue;
                }
		else if(temp[0] == 's')
                {
                        removePrefix(temp);
			tempIndex[i].SA = (int *) malloc(tempIndex[i].length*sizeof(int));
			int j;
			char *number;
			number = strtok(temp," ");
			tempIndex[i].SA[0] = atoi(number);
			for(j = 1; j < tempIndex[i].length; j++)
			{
				number = strtok(NULL," ");
				tempIndex[i].SA[j] = atoi(number);
			}
                }
		else if(temp[0] == 't')
                {
                        removePrefix(temp);
			tempIndex[i].transform = (char *) malloc(strlen(temp)*sizeof(char));
			strcpy(tempIndex[i].transform,temp);
                }
		else if(temp[0] == 'c')
                {
                        removePrefix(temp);
                        int j;
                        char *number;
                        number = strtok(temp," ");
                        tempIndex[i].C[0] = 0;
                        for(j = 1; j < 5; j++)
                        {
                                number = strtok(NULL," ");
                                tempIndex[i].C[j] = atoi(number);
                        }
                }
		else if(temp[0] == 'o')
                {
			int z;
                        char *number;
			removePrefix(temp);
			number = strtok(temp," ");
			tempIndex[i].O[oCount][0] = atoi(number);
			for(z = 1; z < tempIndex[i].length; z++)
			{
                              	number = strtok(NULL," ");
                               	tempIndex[i].O[oCount][z] = atoi(number);
			}
			oCount++;
			if(oCount >= 5)
			{
				oCount = 0;
				i++;
			}
                }
		
	}
	fclose(file);
	return tempIndex;
}

int baseMap(char temp)
{
        if(temp == '$') return 0;
        else if(temp == 'A') return 1;
        else if(temp == 'C') return 2;
        else if(temp == 'G') return 3;
        else if(temp == 'T') return 4;
}

struct output **search(struct query input,int qsc,struct index *interval,int isc)
{
	int z,j;
	struct output **temp = (struct output **) malloc(qsc*sizeof(struct output *));
	for(z = 0; z < qsc; z++)
	{
		temp[z] = (struct output *) malloc(isc*sizeof(struct output));
		for(j = 0; j < isc; j++)
		{
			int i = input.length[z];
			char c = input.sequence[z][0];
			int k = interval[j].C[baseMap(c)] + 1;
			int l = interval[j].C[baseMap(c)+1]; 
			int tempk = k -1;
			while(k <= l && 1 <= i)
			{
				c = input.sequence[z][i-1];
				int x = baseMap(c);
				k = interval[j].C[x] + interval[j].O[x][tempk] + 1;
				l = interval[j].C[x] + interval[j].O[x][l];
				i = i - 1;
			}
			temp[z][j].start = interval[z].SA[k];
			temp[z][j].end = interval[z].SA[l];
			//printf("k= %d\tl= %d\n%d\t%d\n",k,l,interval[0].SA[k],interval[0].SA[l]);
		}
	}	
	return temp;
}

void outputToFile(struct output **out, int qsc, int isc)
{
	FILE *f;
        if(strlen(OUTPUT_FILE) == 0)
        {
                f = fopen("out.bed","w");
        }
        else
        {
                strcat(OUTPUT_FILE,".bed");
                f = fopen(OUTPUT_FILE,"w");
        }
        /*write each instance of M to a file*/
        if (f == NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	int i,j;
	for(i = 0; i < qsc; i++)
	{
		for(j = 0; j < isc; j++)
		{
//			fprintf(f,"%s\t%d\t%d\n",out[i][j].sequence,out[i][j].start,out[i][j].end);
			fprintf(f,"Query: %d\t%d\t%d\n",i+1,out[i][j].start,out[i][j].end);	
		}
	}	
	fclose(f);
}

/*MAIN FUNCTION*/

int main(int argc, char *argv[])
{
	int QseqCount = 0;
	int IseqCount = 0;
	struct query input = manageInputs(argv,argc,&QseqCount);
        struct index *interval = getIndex(&IseqCount);
	//Search
	struct output **out = search(input,QseqCount,interval,IseqCount);
	//Output
	outputToFile(out,QseqCount,IseqCount);
}
