#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "bwp-search.h"
#include <stdbool.h>

bool silent = false;
bool verbose = false;
int MAX_MISMATCHES = 0;
int MAX_LINE_LENGTH = 100000000;
char INTERVAL_FILE[] ="index.bwp";
char OUTPUT_FILE[100];
int GAP = 11;
int EXTENSION = 1;
int subMat[20][20] = {4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2,0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3,-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2,-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3,0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3,-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2,-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2,-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2,-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3,-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2,1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2,0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2,0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2,-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7};


/* MEMORY ALLOCATION FUNCTIONS */

struct input initializeInputStruct(int seqCount, int *seqLength)
{
        struct input query;
        int i;
	int n = 100;
        query.length = (int *) malloc(seqCount * sizeof(int));
        query.name = (char **) malloc(seqCount * sizeof(char *));
        query.sequence = (char **) malloc(seqCount * sizeof(char *));
        for(i = 0; i < seqCount;i++)
        {
                query.name[i] = (char *) malloc(n*(sizeof(char)));
                query.sequence[i] = (char *) malloc(seqLength[i] * sizeof(char));
        }
        return query;
}

/* HELPER FUNCTIONS */

/*  seqCount
 *
 *      Purpose: Find the number of sequences in the query file
 *              Returns: integer
 *              */
int getSeqCount(char *fileName)    /*reutrns number of sequences present in a multi-fasta file*/
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
int *getSeqLength(char *fileName,int seqCount,int charCount) /*returns an array of sequence lengths*/
{
        FILE *file = fopen(fileName,"r");
        char *seq = (char *) malloc(charCount * sizeof(char));
        int *seqLength = (int *) malloc(seqCount * sizeof(int));;
        int i = 0;
        while(fgets(seq,charCount-1,file) != NULL)
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
    	/*m = newArr(seqC,seqLength);*/
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
        if(temp == 0){return 'A';}
        else if(temp == 1){return 'C';}
        else if(temp == 2){return 'D';}
        else if(temp == 3){return 'E';}
	else if(temp == 4){return 'F';}
        else if(temp == 5){return 'G';}
        else if(temp == 6){return 'H';}
	else if(temp == 7){return 'I';}
        else if(temp == 8){return 'K';}
        else if(temp == 9){return 'L';}
	else if(temp == 10){return 'M';}
        else if(temp == 11){return 'N';}
        else if(temp == 12){return 'P';}
	else if(temp == 13){return 'Q';}
        else if(temp == 14){return 'R';}
        else if(temp == 15){return 'S';}
	else if(temp == 16){return 'T';}
        else if(temp == 17){return 'V';}
        else if(temp == 18){return 'W';}
	else if(temp == 19){return 'Y';}
}

int baseMap(char temp)
{
	if(temp == 'A') return 0;
        else if(temp == 'C') return 1;
        else if(temp == 'D') return 2;
        else if(temp == 'E') return 3;
        else if(temp == 'F') return 4;
        else if(temp == 'G') return 5;
        else if(temp == 'H') return 6;
        else if(temp == 'I') return 7;
        else if(temp == 'K') return 8;
        else if(temp == 'L') return 9;
        else if(temp == 'M') return 10;
        else if(temp == 'N') return 11;
        else if(temp == 'P') return 12;
        else if(temp == 'Q') return 13;
        else if(temp == 'R') return 14;
        else if(temp == 'S') return 15;
        else if(temp == 'T') return 16;
        else if(temp == 'V') return 17;
        else if(temp == 'W') return 18;
        else if(temp == 'Y') return 19;
}

void handleS(struct input *query, char *sequence)
{
    	int *seqL = (int *) malloc(sizeof(int));
    	seqL[0] = strlen(sequence);
    	int seqC = 1;
    	*query = initializeInputStruct(seqC,seqL);
    	query->length[0] = seqL[0];
    	strcpy(query->sequence[0],sequence);
	strcpy(query->name[0],"String Input");
}

void handleF(struct input *query,char *fileName)
{
    	int seqCount = getSeqCount(fileName);
    	int charCount = MAX_LINE_LENGTH;
    	int *seqLength = getSeqLength(fileName,seqCount,charCount);
    	*query = initializeInputStruct(seqCount,seqLength);
    	read_fasta(fileName, query);
}

struct input manageInputs(char *argv[], int argc, int *seqCount)
{
	struct input query;
        int c;
        if(argc <= 1)   /*no arguements supplied*/
        {
                printf("Please provide the necessary options\n\nuse -h for usage statement\n");
                exit(0);
        }
        opterr = 0;
	while ((c = getopt (argc, argv, "vShf:s:m:o:d:i:")) != -1) /*options must be added here to be recognized, options followed by : take in a parameter*/
        {
                switch (c)
                {
                        case 'h':
                                printf("\nBurrows Wheeler Nucleotide Alligner\n\nUsage: \"bwp-search <options>\"\n\nOptions:\n\n-f\t\tFor input of a fasta file as a query\n-s\t\tFor input of a string as a query\n-h\t\tFor this usage statement\n-S\t\tTo suppress all output to screen\n-m\t\tTo designate maximum sequence length according to character count\n-o\t\tSpecify output file name (exclude file extentions)\n-d\t\tTo designate the number of allowed mismatches\n-i\t\tTo specify a custom index file\n-v\t\tTo supply output to screen\n\n");
                                exit(0);

                        case 'f':
				if(fileExists(optarg))
                                {
                                	handleF(&query,optarg);
                                	*seqCount = getSeqCount(optarg);
                                	break;
				}
				else
				{
					printf("%s not found, exiting\n",optarg);
                                        exit(0);
				}

                        case 's' :
                                handleS(&query,optarg);
                                *seqCount = 1;
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
				break;
			case 'o' :
                                strcpy(OUTPUT_FILE,optarg);
				break;
			case 'v' :
				verbose = true;
				break;
			case 'd' : 
				MAX_MISMATCHES = atoi(optarg);
				break;
			case 'i' :
				if(fileExists(optarg))
				{	
					strcpy(INTERVAL_FILE,optarg);
				}
				else
				{
					printf("%s is not a valid file\nExiting\n",optarg);
					exit(0);
				}
				break;
			case 'S' :
				silent = true;
				break;
                }
        }
	if(optind < argc)
	{
   		printf ("Non-option arguments supplied\n\nuse -h for usage statement\n");
		exit(0);
	}
        return query;
}

int fileExists(char *temp)
{
	
	FILE *file = fopen(temp, "r");
	if (file)
	{
		fclose(file);
		return 1;
	}
	else
	{
		return 0;
	}
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
		tempIndex[i].O = (int **) malloc(20 * sizeof(int *));
		tempIndex[i].R = (int **) malloc(20 * sizeof(int *));
		for(j = 0; j < 20; j++)
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
		else if(temp[0] == 'q')
		{
			removePrefix(temp);
			tempIndex[i].sequence = (char *) malloc(strlen(temp)*sizeof(char));
			strcpy(tempIndex[i].sequence,temp);
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
		else if(temp[0] == 'f')
                {
                        removePrefix(temp);
                        tempIndex[i].reverse = (char *) malloc(strlen(temp)*sizeof(char));
                        strcpy(tempIndex[i].reverse,temp);
                }
		else if(temp[0] == 'c')
                {
                        removePrefix(temp);
                        int j;
                        char *number;
                        number = strtok(temp," ");
                        tempIndex[i].C[0] = 0;
                        for(j = 1; j < 20; j++)
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
			if(oCount >= 20)
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
                        if(rCount >= 20)
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


int ***calculateD(struct FMidx *index,int isc, struct input query,int qsc)
{
	int i,j,z;
	int ***temp = (int ***) malloc(qsc*sizeof(int **));
	for(i = 0; i < qsc; i++)
	{
		temp[i] = (int **) malloc(isc*sizeof(int *));
		for(z = 0; z < isc; z++)
		{
			temp[i][z] = (int *) malloc(query.length[i]*sizeof(int));
			int low = 1;
			int high = index[z].length - 1;
			int d = 0;
			for(j = 0; j < query.length[i]; j++)
			{
				int k = index[z].C[baseMap(query.sequence[i][j])] + index[z].R[baseMap(query.sequence[i][j])][low-1] + 1;
				int l = index[z].C[baseMap(query.sequence[i][j])] + index[z].R[baseMap(query.sequence[i][j])][high];
				if(index[z].reverse[0] == query.sequence[i][j])
				{
					k = k - 1;
				}
				if(k > l)
				{
					low = 1;
					high = index[z].length - 1;
					d = d + 1;
				}
				temp[i][z][j] = d;
			}
		}
		
	}
	return temp;
}

struct results **inexactSearch(struct input query,int qsc,struct FMidx *index,int isc,int ***D)
{
	int i,j;
	struct results **temp = (struct results **) malloc(qsc*sizeof(struct results *));
	for(i = 0; i < qsc; i++)
	{
		temp[i] = (struct results *) malloc(isc*sizeof(struct results));
		for(j = 0; j < isc; j++)
		{
			temp[i][j].match = NULL;
			char *traceBack = (char *) malloc(sizeof(char) * (MAX_MISMATCHES+query.length[i]));
			temp[i][j].match = inexRecur(index[j],D[i][j],query.sequence[i],query.length[i]-1,MAX_MISMATCHES,1,index[j].length-1,0,1,traceBack,0);
			temp[i][j].match = sortMatches(temp[i][j].match);
		}
	}
	return temp;
}

struct matches *inexRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high,int score,int pState,char *traceBack,int tbIdx)
{
	int tempP = pState;
	int tempS = score;
	struct matches *match = (struct matches *) malloc(sizeof(struct matches));
	match->next = NULL;
	struct matches *results = NULL;
	if(i < 0)
	{
		if(d < D[0]){ return NULL;}	
	}
	else
	{
		if(d < D[i]){ return NULL;}
	}
	if(i < 0)								  
	{
		match->low = low;
		match->high = high;
		match->score = score;
		match->tb = (char *) malloc(sizeof(char) * strlen(traceBack));
		strcpy(match->tb,traceBack);
		match->traceLength = tbIdx;
		return match;
	}
	traceBack[tbIdx] = 'X';
	match = inexRecur(index,D,W,i-1,d-1,low,high,getScore(0,0,tempS,tempP,2),2,traceBack,tbIdx+1); //GAP in index seq
	if(match != NULL){ results = getUnion(match,results);}
	match = NULL;
	tempP = pState;
	tempS = score;
	int j;
	for(j = 0;j < 20; j++)
	{
		int k = index.C[j] + index.O[j][low-1] + 1;
		int l = index.C[j] + index.O[j][high];
		if(index.transform[0] == revBaseMap(j) && low == 1 && (high == index.length-1))
		{
			k = k - 1;
		}
		if(k <= l)
		{
			traceBack[tbIdx] = 'Y';
		        match = inexRecur(index,D,W,i,d-1,k,l,getScore(j,baseMap(W[i]),tempS,tempP,3),3,traceBack,tbIdx+1); //GAP in query
			if(match != NULL){results = getUnion(match,results);}
			match = NULL;
			tempP = pState;
			tempS = score;
			if(revBaseMap(j) == W[i]) //match
			{
				traceBack[tbIdx] = 'M';
				match = inexRecur(index,D,W,i-1,d,k,l,getScore(j,baseMap(W[i]),tempS,tempP,1),1,traceBack,tbIdx+1);
				if(match != NULL){ results = getUnion(match,results);}
				match = NULL;
				tempP = pState;
				tempS = score;
			}
			else //mismatch
			{
				traceBack[tbIdx] = 'U';
				match = inexRecur(index,D,W,i-1,d-1,k,l,getScore(j,baseMap(W[i]),tempS,tempP,1),1,traceBack,tbIdx+1);
				if(match != NULL){ results = getUnion(match,results);}
				match = NULL;
				tempP = pState;
				tempS = score;
			}
		}
	}
	return results;
}

struct matches *getUnion(struct matches *head1, struct matches *head2)
{
    	struct matches *result = NULL;
    	struct matches *t1 = head1, *t2 = head2;
	while (t1 != NULL)
    	{
        	push(&result, t1->low,t1->high,t1->score,t1->tb,t1->traceLength); //put new matches at top of stack
        	t1 = t1->next;
    	}
	while (t2 != NULL)
    	{
//      		if(!isPresent(result, t2->low,t2->high))
//		{
            		push(&result, t2->low,t2->high,t2->score,t2->tb,t2->traceLength);
//		}
//		else //present but checking if score is better
//		{
//			getHighScore(result,t2->low,t2->high,t2->score,t2->tb,t2->traceLength);
//		}
        	t2 = t2->next;
    	}
    	return result;
}

void getHighScore(struct matches *head, int k, int l, int score,char *traceback,int traceLength)
{
	struct matches *t = head;
    	while (t != NULL)
    	{
        	if((t->low == k)&&(t->high == l)&&(t->score < score))
        	{
                	t->score = score;
			strcpy(t->tb,traceback);
			t->traceLength = traceLength;
			return;
        	}
        	t = t->next;
    	}	
}

void push(struct matches** head_ref, int k,int l,int score,char *traceback,int traceLength)
{
    	struct matches *new_node = (struct matches *) malloc(sizeof(struct matches));
    	new_node->low = k;
	new_node->high = l;
	new_node->score = score;
	new_node->tb = (char *) malloc(sizeof(char) * strlen(traceback));
	strcpy(new_node->tb,traceback);
	new_node->traceLength = traceLength;
    	new_node->next = (*head_ref);
    	(*head_ref) = new_node;
}

int isPresent(struct matches *head, int k, int l)
{
    struct matches *t = head;
    while (t != NULL)
    {
    	if((t->low == k)&&(t->high == l))
	{
        	return 0; //should be 1 
	}
        t = t->next;
    }
    return 0;
}

struct matches *pointToTail(struct matches *match)
{
	if(match != NULL)
        {
                while(match->next != NULL)
                {
                        match = match->next;
                }
        }
	return match;
}

struct matches *sortMatches(struct matches *match)
{
	mergeSortMatches(&match);
	//Potentially more here
	return match;	
}

void mergeSortMatches(struct matches** headRef)
{
    struct matches* head = *headRef;
    struct matches* a;
    struct matches* b;
    if ((head == NULL) || (head -> next == NULL))
    {
        return;
    }
    frontbacksplit(head, &a, &b);
    mergeSortMatches(&a);
    mergeSortMatches(&b);
    *headRef = sortedmergeMatch(a, b);
}
 
struct matches* sortedmergeMatch(struct matches* a, struct matches* b)
{
    struct matches* result = NULL;
 
    if (a == NULL)
        return(b);
    else if (b == NULL)
        return(a);
 
    if ( a->score >= b->score)
    {
        result = a;
        result->next = sortedmergeMatch(a->next, b);
    }
    else
    {
        result = b;
        result->next = sortedmergeMatch(a, b->next);
    }
    return(result);
}
 
void frontbacksplit(struct matches* source, struct matches** frontRef, struct matches** backRef)
{
    struct matches* fast;
    struct matches* slow;
    if (source==NULL || source->next==NULL)
    {
        *frontRef = source;
        *backRef = NULL;
    }
    else
    {
        slow = source;
        fast = source -> next;
        while (fast != NULL)
        {
            fast = fast -> next;
            if (fast != NULL)
            {
                slow = slow -> next;
                fast = fast -> next;
            }
	}
    }
 
    *frontRef = source;
    *backRef = slow -> next;
    slow -> next = NULL;
}

int getScore(int l1, int l2,int score, int p, int c)
{
	int temp = score;
	if(p <= 1)
	{
		if(c == 1) //match
                {
                        temp = temp + subMat[l1][l2];
                        return temp;
                }
                else if(c == 2) //GAP in index
                {
                        temp = temp - GAP;
                        return temp;
                }
                else if(c == 3) //GAP in query
                {
                        temp = temp - GAP;
			return temp;
                }
	}
	else if(p == 2)
	{
		if(c == 1)
                {
                        temp = temp + subMat[l1][l2];
                        return temp;
                }
                else if(c == 2) //extension
                {
                        temp = temp - EXTENSION;
                        return temp;
                }
                else if(c == 3) 
                {
                        temp = temp - GAP;
                        return temp;
                }
	}
	else if(p == 3)
	{
		if(c == 1)
                {
                        temp = temp + subMat[l1][l2];
                        return temp;
                }
                else if(c == 2)
                {
                        temp = temp - GAP;
                        return temp;
                }
                else if(c == 3) //extension
                {
                        temp = temp - EXTENSION;
                        return temp;
                }
	}
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

void outputToFile(struct results **out, int qsc, int isc,struct FMidx *index,struct input query)
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
//      write each instance of M to a file
        if (f == NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
        int i,j;
	fprintf(f,"Database\tStart Pos\tEnd Pos\tQuery\tScore\n----------------------------------------------------------------\n");
        for(i = 0; i < qsc; i++)
        {
                for(j = 0; j < isc; j++)
                {
			struct matches *temp = (struct matches *) malloc(sizeof(struct matches));
                        temp = out[i][j].match;
                        if(temp != NULL)
                        {
                                while(temp != NULL)
                                {
					fprintf(f,"%s\t",index[j].desc);
                                        fprintf(f,"%d\t%d\t%s\t%d\n",index[j].SA[temp->low],query.length[i]+index[j].SA[temp->low]-1,query.name[i],temp->score);
                                        temp = temp->next;
                                }
                        }
                }
        }
        fclose(f);
}

void printInResults(struct results **out,int qsc,int isc,struct FMidx *index, struct input query)
{
	int i,j,z;
	int miss = 0;
	for(i = 0; i < qsc; i++)
	{
		for(j = 0; j < isc; j++)
		{
			struct matches *temp = (struct matches *) malloc(sizeof(struct matches));
			temp = out[i][j].match; //causing seg fault
			printf("\nindex: %s\nquery: %s\n\n-----------------------------------------------\n\n",index[j].desc,query.name[i]);
			if(temp!=NULL)
			{
				while(temp!=NULL)
				{
//					printf("tb:%s\n",temp->tb);
					printf("Score: %d\n\nIndex\t",temp->score);
					for(z = 0; z < temp->traceLength; z++)
					{
						if(z % 59 == 0 && z != 0)
						{
							printf("\n");
						}
						else
						{
							if(temp->tb[z] == 'M' || temp->tb[z] == 'U')
							{
								printf("%c",index[j].sequence[index[j].SA[temp->low]+(z-1-miss)]);
							}
							else if(temp->tb[z] == 'X')
							{
								printf("-");
								miss++;
							}
							else if(temp->tb[z] == 'Y')
                                                	{
								printf("%c",index[j].sequence[index[j].SA[temp->low]+(z-1-miss)]);
                                                	}	
						}
					}
					printf("\n\t");
					int missX = 0;
					int missY = 0;
					for(z = 0; z < temp->traceLength; z++)
                                        {
						if(z % 59 == 0 && z != 0)
                                                {
                                                        printf("\n");
                                                }
                                                else
                                                {
							if(query.sequence[i][z-missY] == index[j].sequence[index[j].SA[temp->low]+z-1-missX] && (temp->tb[z] != 'X' || temp->tb[z] != 'Y'))
							{
								printf("|");
							}
							else if(temp->tb[z] == 'X')
							{
								printf(" ");
								missX++;
							}
							else if( temp->tb[z] == 'Y')
							{
								printf(" ");
                                                                missY++;
							}
							else
							{
								printf(" ");
							}
						}
					}
					printf("\nQuery\t");
					miss = 0;
					for(z = 0; z < temp->traceLength; z++)
                                        {
						if(z % 59 == 0 && z != 0)
                                                {
                                                        printf("\n");
                                                }
                                                else
                                                {
	                                                if(temp->tb[z] == 'M' || temp->tb[z] == 'U')
	                                                {
	                                                        printf("%c",query.sequence[i][z-miss]);
	                                                }
	                                                else if(temp->tb[z] == 'X')
	                                                {
								printf("%c",query.sequence[i][z-miss]);
	                                                }
	                                                else if(temp->tb[z] == 'Y')
	                                                {
								printf("-");
								miss++;
	                                                }
						}
                                        }
					printf("\n\nStart: %d\nEnd: %d\n\n-----------------------------------------------\n\n",index[j].SA[temp->low],index[j].SA[temp->low]+temp->traceLength-1);
					temp = temp->next;
				}
			}
			else
			{
				printf("No Matches Found\n");
			}
		}
	}	
}

void printKLs(struct results **out, int qsc, int isc)
{
	int i,j;
	for(i = 0; i < qsc; i ++)
	{
		for(j = 0; j < isc; j++)
		{
			printf("index: %d\tquery: %d\n",j+1,i+1);
			struct matches *temp = (struct matches *) malloc(sizeof(struct matches));
                        temp = out[i][j].match;
			if(temp!=NULL)
                        {
                                while(temp!=NULL)
                                {
					printf("k: %d\tl: %d\n",temp->low,temp->high);
				}
			}
		}
	}
}

/*MAIN FUNCTION*/

int main(int argc, char *argv[])
{
	int QseqCount = 0;
	int IseqCount = 0;
	struct input query = manageInputs(argv,argc,&QseqCount);
        struct FMidx *index = getIndex(&IseqCount);
	int ***D = calculateD(index,IseqCount,query,QseqCount);
	//Search
//	struct output **out = search(query,QseqCount,index,IseqCount);
	struct results **out = inexactSearch(query,QseqCount,index,IseqCount,D);
	//Output
	outputToFile(out,QseqCount,IseqCount,index,query);
//	printResults(out,QseqCount,IseqCount,index);
	if(verbose)
	{
		printKLs(out,QseqCount,IseqCount);
	}
	if(!silent)
	{
		printInResults(out,QseqCount,IseqCount,index,query);
	}
}
