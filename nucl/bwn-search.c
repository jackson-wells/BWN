#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "bwn-search.h"

int MAX_LINE_LENGTH = 10000000;
#define INTERVAL_FILE "index.bwn"
char OUTPUT_FILE[100];

/* MEMORY ALLOCATION FUNCTIONS */

struct input initializeInputStruct(int sCount, int *sLength)
{
        struct input query;
        int i;
        int n = 100;
        query.length = (int *) malloc(sCount * sizeof(int));
        query.name = (char **) malloc(sCount * sizeof(char *));
        query.sequence = (char **) malloc(sCount * sizeof(char *));
        for(i = 0; i < sCount;i++)
        {
                query.name[i] = (char *) malloc(n*(sizeof(char)));
                query.sequence[i] = (char *) malloc(sLength[i] * sizeof(char));
        }
        return query;
}

/* HELPER FUNCTIONS */

/*  seqCount
 *
 *      Purpose: Find the number of sequences in the query file
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
 *      Purpose: get the lengths of query sequences from file
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
                seqLength[i] = strlen(seq);
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
void read_fasta(char *fileName, struct input *query)
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
                    strcpy(query->name[i],temp);
            }
            else if(temp[0] == '\n') /*if line is empty*/
            {
                    continue;
            }
            else    /*if line contains a nucleotide sequence*/
            {
                    strtok(temp,"\n"); /*strings read from file have extra \n added by file read*/
                    //strcat(temp,"$");
                    strcpy(query->sequence[i],temp);    /*saving string in memory*/
                    query->length[i] = strlen(temp);
                    i++;
            }
    }
    fclose(file);
        /*free(temp);*/
}

char revBaseMap(int temp)
{
        if(temp == 0){ return '$';}
        else if(temp == 1){return 'A';}
        else if(temp == 2){return 'C';}
        else if(temp == 3){return 'G';}
        else if(temp == 4){return 'T';}
}

int baseMap(char temp)
{
        if(temp == '$') return 0;
        else if(temp == 'A') return 1;
        else if(temp == 'C') return 2;
        else if(temp == 'G') return 3;
        else if(temp == 'T') return 4;
}

void handleS(struct input *query, char *sequence)
{
    printf("\nHere is the string you entered:\n%s\n\n",sequence);
    int *seqL = (int *) malloc(sizeof(int));
    seqL[0] = strlen(sequence);
    int seqC = 1;
//    strcat(sequence,"$");
    *query = initializeInputStruct(seqC,seqL);
    query->length[0] = seqL[0];
    strcpy(query->sequence[0],sequence);
    strcpy(query->name[0],"Command Line String Input");
}

void handleF(struct input *query,char *fileName)
{
    	int sCount = seqCount(fileName);
    	int cCount = MAX_LINE_LENGTH;
    	int *sLength = seqLength(fileName,sCount,cCount);
    	*query = initializeInputStruct(sCount,sLength);
    	read_fasta(fileName, query);
}

struct input manageInputs(char *argv[], int argc, int *sCount)
{
	struct input query;
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
                                printf("\nBurrows Wheeler Nucleotide Alligner\n\nUsage: \"bwn-search <options>\"\n\nOptions:\n\n-f\t\tFor query of a fasta file\n-s\t\tFor query of a string\n-h\t\tFor this usage statement\n-m\t\tTo designate maximum sequence length according to character count\n-o\t\tSpecify output file name (exclude file extentions)\n\n");
                                exit(0);

                        case 'f':
                                handleF(&query,optarg);
                                *sCount = seqCount(optarg);
                                break;

                        case 's' :
                                handleS(&query,optarg);
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
        return query;
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

struct FMidx *getIndex(int *IseqCount)
{
	int seqCount = getCount();
	int *seqLength = getLength(seqCount);
	*IseqCount = seqCount;
	int charCount = MAX_LINE_LENGTH;
	char *temp = (char *) malloc((MAX_LINE_LENGTH+1)*sizeof(char));
	struct FMidx *tempIndex = (struct FMidx *) malloc(seqCount * sizeof(struct FMidx));
	FILE *file = fopen(INTERVAL_FILE,"r");
    	int i,j;
	for(i = 0; i < seqCount; i++)
	{
		tempIndex[i].O = (int **) malloc(5 * sizeof(int *));
		tempIndex[i].R = (int **) malloc(5 * sizeof(int *));
		for(j = 0; j < 5; j++)
		{
			tempIndex[i].R[j] = (int *) malloc(seqLength[i] * sizeof(int));
			tempIndex[i].O[j] = (int *) malloc(seqLength[i] * sizeof(int));
		}
		tempIndex[i].length = seqLength[i];
	}
	i = 0;
	int oCount = 0;
	int rCount = 0;
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
                        for(j = 1; j < 6; j++)
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
			}
                }
		else if(temp[0] == 'r')
                {
                        int z;
                        char *number;
                        removePrefix(temp);
                        number = strtok(temp," ");
                        tempIndex[i].R[rCount][0] = atoi(number);
                        for(z = 1; z < tempIndex[i].length; z++)
                        {
                                number = strtok(NULL," ");
                                tempIndex[i].R[rCount][z] = atoi(number);
                        }
                        rCount++;
                        if(rCount >= 5)
                        {
                                rCount = 0;
                                i++;
                        }
                }
		
	}
	fclose(file);
	return tempIndex;
}

struct output **search(struct input query,int qsc,struct FMidx *index,int isc)
{
	int z,j;
	struct output **temp = (struct output **) malloc(qsc*sizeof(struct output *));
	for(z = 0; z < qsc; z++)
	{
		temp[z] = (struct output *) malloc(isc*sizeof(struct output));
		for(j = 0; j < isc; j++)
		{
			temp[z][j].sequence = (char *) malloc(query.length[z]*sizeof(char));
			strcpy(temp[z][j].sequence,query.sequence[z]);
			int i = query.length[z]-1;
			char c = query.sequence[z][i];
			int low = index[j].C[baseMap(c)] + 1;
			int high = index[j].C[baseMap(c)+1];
			while(low <= high && 1 <= i)
			{
				c = query.sequence[z][i-1];
				low = index[j].C[baseMap(c)] + index[j].O[baseMap(c)][low-2] +1; 
				high = index[j].C[baseMap(c)] + index[j].O[baseMap(c)][high-1]; //-1 for 0-base
				i = i - 1;
			}
			if(high < low)
			{
				temp[z][j].low = 0;
				temp[z][j].high = 0;
			}
			else
			{
				temp[z][j].low = low-1;
				temp[z][j].high = high-1;//index[j].SA[high-1];
			}
		}
	}	
	return temp;
}


int **calculateD(struct FMidx *index,int isc, struct input query,int qsc)
{
	int i,j;
	int **temp = (int **) malloc(qsc*sizeof(int *));
	for(i = 0; i < qsc; i++)
	{
		temp[i] = (int *) malloc(query.length[i]*sizeof(int));
		int low = 1;
		int high = index[i].length;
		int d = 0;
		for(j = 0; j < query.length[i]; j++)
		{
			low = index[i].C[baseMap(query.sequence[i][j])] + index[i].R[baseMap(query.sequence[i][j])][low-2] +1;
			high = index[i].C[baseMap(query.sequence[i][j])] + index[i].R[baseMap(query.sequence[i][j])][high-1];
			if(low > high)
			{
				low = 1;
				high = index[i].length;
				d = d + 1;
			}
			temp[i][j] = d;
		}
	}
	return temp;
}

struct output ***inExactSearch(struct input query,int qsc,struct FMidx *index,int isc,int **D)
{
	int i,j;
	struct output ***temp = (struct output ***) malloc(qsc*sizeof(struct output **));
	for(i = 0; i < qsc; i++)
	{
		temp[i] = (struct output **) malloc(isc*sizeof(struct output *));
		for(j = 0; j < isc; j++)
		{
			temp[i][j] = inexRecur(index[i],D[i],query.sequence[i],query.length[i],D[i][j],1,index[i].length);
		}
	}
}

struct output *inexRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high)
{
	if(d < D[i]){ return NULL;}
	struct output *temp =  (struct output *) malloc(5*sizeof(struct output)); //5 is chosen # of allowed 
	if(i < 1)								  //matches 
	{ 
		temp[i].low = low;
		temp[i].high = high;
		return temp;
	}
	temp = inexRecur(index,D,W,i-1,d-1,low,high); //overwrite temp?
	int j;
	for(j = 1;j < 5; j++)
	{
		low = index.C[j] + index.O[j][low-2] + 1;
		high = index.C[j] + index.O[j][high-1];
		if(low <= high)
		{
			temp = inexRecur(index,D,W,i,d-1,low,high);
			if(revBaseMap(j) == W[i])
			{
				//for(i = 0; i < 5; i++)
				//{
				//	temp2 = temp;
				//	
				//	compare returned valuses? 	
				//}
				temp = inexRecur(index,D,W,i-1,d,low,high);
			}
			else
			{
				temp = inexRecur(index,D,W,i-1,d-1,low,high);
			}
		}
	}
	return temp;
}

void printResults(struct output **out, int qsc, int isc,struct FMidx *index)
{
	int i,j;
	for(i = 0; i < qsc; i++)
	{
		printf("Query Sequence: %d\n",i+1);
		for(j = 0; j < isc; j++)
		{
			printf("\tDatabase Sequence: %d\n",j+1);
			if(out[i][j].low == 0 && out[i][j].high == 0)
			{
				printf("\t\tNot Found\n\n");
			}
			else
			{
				if(out[i][j].high == out[i][j].low )
				{
					printf("\t\tFound \"%s\" starting at the %dth position\n\n",out[i][j].sequence,index[j].SA[out[i][j].low]);
				}
				else
				{
					int temp = (out[i][j].high - out[i][j].low)+1;
					int z;
					printf("\t\tFound \"%s\" starting at positions:\n",out[i][j].sequence);
					for(z = 0; z < temp; z++)
					{
						printf("\t\t%d\n",index[j].SA[out[i][j].low+z]);	
					}
					printf("\n");
				}
			}
		}
		printf("\n");
	}	
}

void outputToFile(struct output **out, int qsc, int isc,struct FMidx *index,struct input query)
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
	fprintf(f,"Query\tDatabase\tStart Position\n--------------------------------------\n");
        for(i = 0; i < qsc; i++)
        {
                fprintf(f,"%s\t",query.name[i]);
                for(j = 0; j < isc; j++)
                {
			if(j > 0){ fprintf(f,"\t");}
                        fprintf(f,"%s\t\t",index[j].desc);
                        if(out[i][j].low == 0 && out[i][j].high == 0)
                        {
                                fprintf(f,"Not Found\n");
                        }
                        else
                        {
                                if(out[i][j].high == out[i][j].low)
                                {
                                        fprintf(f,"%d\n",index[j].SA[out[i][j].low]);
                                }
                                else
                                {
                                        int temp = (out[i][j].high - out[i][j].low)+1;
                                        int z;
                                        for(z = 0; z < temp; z++)
                                        {
                                                fprintf(f,"%d",index[j].SA[out[i][j].low+z]);
						if(z < (temp-1)){ fprintf(f,","); }
                                        }
                                        fprintf(f,"\n");
                                }
                        }
                }
                fprintf(f,"\n");
        }
        fclose(f);
}


/*MAIN FUNCTION*/

int main(int argc, char *argv[])
{
	int QseqCount = 0;
	int IseqCount = 0;
	struct input query = manageInputs(argv,argc,&QseqCount);
        struct FMidx *index = getIndex(&IseqCount);
	//int **D = calculateD(index,IseqCount,query,QseqCount);
	//Search
	struct output **out = search(query,QseqCount,index,IseqCount);
	//struct output ***out = inExactSearch(query,QseqCount,index,IseqCount,D);
	//Output
	//outputToFile(out,QseqCount,IseqCount,index,query);
	printResults(out,QseqCount,IseqCount,index);
}
