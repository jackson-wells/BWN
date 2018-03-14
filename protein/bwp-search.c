#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "bwp-search.h"
#include <stdbool.h>
#include <getopt.h>

bool silent = false;
bool verbose = false;
int MAX_MISMATCHES = 0;
int MAX_LINE_LENGTH = 100000000;
char INTERVAL_FILE[] ="index.bwp";
char OUTPUT_FILE[100];
int GAP = 11;
int EXTENSION = 1;
int subMat[20][20] = {{4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2},{0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2},{-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3},{-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2},{-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3},{0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3},{-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2},{-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1},{-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2},{-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1},{-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1},{-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2},{-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3},{-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1},{-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2},{1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2},{0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2},{0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1},{-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2},{-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7}};


/* 					*
 * 	MEMORY ALLOCATION FUNCTIONS 	*
 * 					*/



/*
	initializeInputStruct
	
	allocates memory for the variables contained in the input
	structure object

	seqCount: integer variable containing the number of seqeunces presint
                in the input file

        seqLength: an array of integers containing the lengths of each sequence
                in the input file

	returns: a structure variable prepared to contain information from the
		input file
 
*/
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

/* 				*
 * 	HELPER FUNCTIONS 	*
 * 				*/


/*  
	seqCount

      	computes the number of sequences in the query file

	fileName: string variable containing the name of the input file	

        returns: an integer variable containing the number of 
		sequence in the input file

*/
int getSeqCount(char *fileName) 
{
        FILE *file = fopen(fileName,"r");
        int lineCount = 0;
        char *temp = (char *) malloc((MAX_LINE_LENGTH+1) * sizeof(char));
        while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)
        {
		if(strstr(temp,">"))
		{
                	lineCount++;
		}
        }
        fclose(file);
        free(temp);
        return lineCount;
}


/* 
	removePrefix

 	removes the first two characters from a string

	input: character array containing a sequence with a prefix

 	returns: character array without its first two characters

*/
char *removePrefix(char *input)
{
	memmove(input, input+2, strlen(input));
	return input;
}


/*  
	getSeqLength

      	get the lengths of query sequences from input file

	fileName: string variable containing the name of the input file

	seqCount: integer variable containing the number of query sequences
              
	returns: array of integers

*/
int *getSeqLength(char *fileName,int seqCount) /*returns an array of sequence lengths*/
{
	int i = 0;
        FILE *file = fopen(fileName,"r");
        char *seq = (char *) malloc(MAX_LINE_LENGTH * sizeof(char));
        int *seqLength = (int *) malloc(seqCount * sizeof(int));;
        while(fgets(seq,MAX_LINE_LENGTH-1,file) != NULL)
        {
            	if(seq[0] != '>' && seq[0] != '\n')
            	{
            	    	seqLength[i] = strlen(seq);
            	    	i++;
            	}
        }
        fclose(file);
        free(seq);
        return seqLength;
}


/*
 
	getLength

	retrieves the lengths of the index sequences

	seqCount: integer varaible containing the number of index sequences

	returns: array of integers

*/
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
	free(temp);
        return length;
}

/*  
	read_fasta

      	stores FASTA file information into memory

	fileName: string variable containing the name of the input file

	query: array of structures to contain input sequence information
	
*/
void read_fasta(char *fileName, struct input *query)
{
    	char *temp = (char *) malloc(MAX_LINE_LENGTH+1 * sizeof(char));
    	/*m = newArr(seqC,seqLength);*/
    	FILE *file = fopen(fileName,"r");
    	int i = 0;
    	while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)   /*loops through each line of a file*/
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
                    	/*strcat(temp,"$");*/
                    	strcpy(query->sequence[i],temp);    /*saving string in memory*/
                    	query->length[i] = strlen(temp);
                    	i++;
            	}
    	}
    	fclose(file);
        free(temp);
}


/*
 
	reverse

	computes the reverse of a supplied string

	str: character array to be reversed

	returns: a reversed version of the input string

*/
char *reverse(char *str)
{
	int i,len,end;
	char *tStr;
	len = strlen(str);
	end = len-1;
	tStr = (char *) malloc(sizeof(char) * len);
	for( i=0 ; i< len ; i++)
	{
		tStr[i] = str[end];
		end--;
	}
	return tStr;
}

/*
 
	revBaseMap

	handles input of integers that corressponds to a character 
	from the amino acid alphabet

        temp: integer  variable containing a value that corresponds to
		a character from the amino acid alphabet

	returns: character corresponding to the supplied integer

*/
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
	else
        {
                printf("Invalid integer '%d'\n",temp);
                exit(1);
        }
}


/*
 
	baseMap

	handles input of characters from the amino acid alphabet

        temp: character variable containing a capitol letter from
                the amino acid alphabet

	returns: 
		if in alphabet: corresponding integer for the input character

                if not in alphabet: error message containing the invalid
                character and exits with status 1

*/
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
	else
        {
                printf("Invalid character '%c'\n",temp);
                exit(1);
        }
}

/*
 
	fileExists
	
	verify that a supplied file exists and is openable

	temp: character array containing the name of a file

	retunrs:
		0: if the file does not exists
	
		1: if the file does exist

*/
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
                return 0; /*should be 1 */
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

/*
 
	getCount

	gets the number of sequences contained in the index file
	
	returns: an integer variable corresponding to the number
		of sequences inside the index file

*/
int getCount(void)
{
	char temp[10];
	FILE *file = fopen(INTERVAL_FILE,"r");
	fgets(temp,10,file);
	fclose(file);
	removePrefix(temp);
	return atoi(temp);
}

struct matches *getUnion(struct matches *head1, struct matches *head2)
{
        struct matches *result = NULL;
        struct matches *t1 = head1, *t2 = head2;
        while (t1 != NULL)
        {
                push(&result, t1->low,t1->high,t1->score,t1->tb,t1->traceLength); /*put new matches at top of stack*/
                t1 = t1->next;
        }
        while (t2 != NULL)
        {
/*              if(!isPresent(result, t2->low,t2->high))
 *              {*/
                push(&result, t2->low,t2->high,t2->score,t2->tb,t2->traceLength);
/*              }*/
                t2 = t2->next;
        }
        return result;
}

/*					*
 *	INPUT HANDLING FUNCTIONS	*
 *					*/

/*

      	manageInputs

        handles command-line input from user

        argc: integer variable containing the number of arguements passed in

        argv: 2d character array containing the arguement strings

        seqCount: integer array to hold the number of sequences to be processed
        	in the run of the program

	returns: a structure variable containing input information
                                                                 */
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

/*

        handleS

       	handles sequence information supplied via command-line

       	query: structure variable containing input information

	sequence: protein sequence input by the user

*/
void handleS(struct input *query, char *sequence)
{
        int seqC = 1;
        int *seqL = (int *) malloc(sizeof(int));
        seqL[0] = strlen(sequence);
        *query = initializeInputStruct(seqC,seqL);
        query->length[0] = seqL[0];
        strcpy(query->sequence[0],sequence);
        strcpy(query->name[0],"String Input");
}


/*

       	handleF

     	handles sequence information supplied via files

        query: structure variable containing file information

	fileName: character string containing the name of the input file

*/
void handleF(struct input *query,char *fileName)
{
        int seqCount = getSeqCount(fileName);
        int *seqLength = getSeqLength(fileName,seqCount);
        *query = initializeInputStruct(seqCount,seqLength);
        read_fasta(fileName, query);
}

/*
 
	getIndex



	IseqCount:

	returns:

*/
struct FMidx *getIndex(void)
{
	int j,z,i,rCount,oCount;
	char *number;
	char *temp;
	FILE *file;
	struct FMidx *tempIndex;
	int seqCount = getCount();
	int *seqLength = getLength(seqCount);
	temp = (char *) malloc((MAX_LINE_LENGTH+1)*sizeof(char));
	tempIndex = (struct FMidx *) malloc(seqCount * sizeof(struct FMidx));
	file = fopen(INTERVAL_FILE,"r");
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
	oCount = 0;
	rCount = 0;
    	while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)   /*loops through each line of a file*/
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

/*				*
 *	SORTING FUCNTIONS	*
 *				*/

struct matches *sortMatches(struct matches *match)
{
        mergeSortMatches(&match);
        /*Potentially more here*/
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

/*				*
 *	SEARCH FUNCTIONS	*
 *				*/

struct output **search(struct input query,int qsc,struct FMidx *index,int isc)
{
	int z,j,i,low,high;
	char c;
	struct output **temp = (struct output **) malloc(qsc*sizeof(struct output *));
	for(z = 0; z < qsc; z++)
	{
		temp[z] = (struct output *) malloc(isc*sizeof(struct output));
		for(j = 0; j < isc; j++)
		{
			temp[z][j].sequence = (char *) malloc(query.length[z]*sizeof(char));
			strcpy(temp[z][j].sequence,query.sequence[z]);
			i = query.length[z]-1;
			c = query.sequence[z][i];
			low = index[j].C[baseMap(c)] + 1;
			high = index[j].C[baseMap(c)+1];
			while(low <= high && 1 <= i)
			{
				c = query.sequence[z][i-1];
				low = index[j].C[baseMap(c)] + index[j].O[baseMap(c)][low-2] +1; 
				high = index[j].C[baseMap(c)] + index[j].O[baseMap(c)][high-1]; /*-1 for 0-base*/
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
				temp[z][j].high = high-1;/*index[j].SA[high-1];*/
			}
		}
	}	
	return temp;
}


int ***calculateD(struct FMidx *index,int isc, struct input query,int qsc)
{
	int i,j,z;
	int ***temp = (int ***) malloc(qsc*sizeof(int **));
	int low = 1;
	int d = 0;
	int high;
	int k,l;
	for(i = 0; i < qsc; i++)
	{
		temp[i] = (int **) malloc(isc*sizeof(int *));
		for(z = 0; z < isc; z++)
		{
			temp[i][z] = (int *) malloc(query.length[i]*sizeof(int));
			low = 1;
			high = index[z].length - 1;
			d = 0;
			for(j = 0; j < query.length[i]; j++)
			{
				k = index[z].C[baseMap(query.sequence[i][j])] + index[z].R[baseMap(query.sequence[i][j])][low-1] + 1;
				l = index[z].C[baseMap(query.sequence[i][j])] + index[z].R[baseMap(query.sequence[i][j])][high];
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
	char *traceBack;
	for(i = 0; i < qsc; i++)
	{
		temp[i] = (struct results *) malloc(isc*sizeof(struct results));
		for(j = 0; j < isc; j++)
		{
			temp[i][j].match = NULL;
			traceBack = (char *) malloc(sizeof(char) * (MAX_MISMATCHES+query.length[i]));
			temp[i][j].match = inexRecur(index[j],D[i][j],query.sequence[i],query.length[i]-1,MAX_MISMATCHES,1,index[j].length-1,0,1,traceBack,0);
			temp[i][j].match = sortMatches(temp[i][j].match);
		}
	}
	return temp;
}

struct matches *inexRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high,int score,int pState,char *traceBack,int tbIdx)
{
	int j;
	char *tempTB = (char *) malloc(sizeof(char)*tbIdx);
	int tempP = pState;
        int tempS = score;
        struct matches *match = (struct matches *) malloc(sizeof(struct matches));
	struct matches *results = NULL;
	strcpy(tempTB,traceBack);
	match->next = NULL;
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
		match->tb = (char *) malloc(sizeof(char) * tbIdx);
		strncpy(match->tb,tempTB,tbIdx);
		strcpy(match->tb,reverse(match->tb));
		match->traceLength = strlen(match->tb);
		return match;
	}
	tempTB[tbIdx] = 'X';
	match = inexRecur(index,D,W,i-1,d-1,low,high,getScore(0,0,tempS,tempP,2),2,tempTB,tbIdx+1); /*GAP in index seq*/
	if(match != NULL){ results = getUnion(match,results);}
	match = NULL;
	tempP = pState;
	tempS = score;
	for(j = 0; j < 20; j++)
	{
		int k = index.C[j] + index.O[j][low-1] + 1;
		int l = index.C[j] + index.O[j][high];
		if(index.transform[0] == revBaseMap(j) && low == 1 && (high == index.length-1))
		{
			k = k - 1;
		}
		if(k <= l)
		{
			tempTB[tbIdx] = 'Y';
		        match = inexRecur(index,D,W,i,d-1,k,l,getScore(j,baseMap(W[i]),tempS,tempP,3),3,tempTB,tbIdx+1); /*GAP in query*/
			if(match != NULL){results = getUnion(match,results);}
			match = NULL;
			tempP = pState;
			tempS = score;
			if(revBaseMap(j) == W[i]) /*match*/
			{
				tempTB[tbIdx] = 'M';
				match = inexRecur(index,D,W,i-1,d,k,l,getScore(j,baseMap(W[i]),tempS,tempP,1),1,tempTB,tbIdx+1);
				if(match != NULL){ results = getUnion(match,results);}
				match = NULL;
				tempP = pState;
				tempS = score;
			}
			else /*mismatch*/
			{
				tempTB[tbIdx] = 'U';
				if(subMat[j][baseMap(W[i])] > 0)
				{
					match = inexRecur(index,D,W,i-1,d,k,l,getScore(j,baseMap(W[i]),tempS,tempP,1),1,tempTB,tbIdx+1);
				}
				else
				{
					match = inexRecur(index,D,W,i-1,d-1,k,l,getScore(j,baseMap(W[i]),tempS,tempP,1),1,tempTB,tbIdx+1);
				}
				if(match != NULL){ results = getUnion(match,results);}
				match = NULL;
				tempP = pState;
				tempS = score;
			}
		}
	}
	return results;
}

int getScore(int l1, int l2,int score, int p, int c)
{
	int temp = score;
	if(p <= 1)
	{
		if(c == 1) /*match*/
                {
                        temp = temp + subMat[l1][l2];
                        return temp;
                }
                else if(c == 2) /*GAP in index*/
                {
                        temp = temp - GAP;
                        return temp;
                }
                else if(c == 3) /*GAP in query*/
                {
                        temp = temp - GAP;
			return temp;
                }
		else
                {
                        printf("Invalid current state value '%d'\n",c);
                        exit(0);
                }
	}
	else if(p == 2)
	{
		if(c == 1)
                {
                        temp = temp + subMat[l1][l2];
                        return temp;
                }
                else if(c == 2) /*extension*/
                {
                        temp = temp - EXTENSION;
                        return temp;
                }
                else if(c == 3) 
                {
                        temp = temp - GAP;
                        return temp;
                }
		else
                {
                        printf("Invalid current state value '%d'\n",c);
                        exit(0);
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
                else if(c == 3) /*extension*/
                {
                        temp = temp - EXTENSION;
                        return temp;
                }
		else
        	{
        	        printf("Invalid current state value '%d'\n",c);
                	exit(0);
        	}
	}
	else
	{
		printf("Invalid previous state value '%d'\n",p);
		exit(0);
	}
}

/*				*
 *	OUTPUT FUNCTIONS	*
 *				*/

void outputToFile(struct results **out, int qsc, int isc,struct FMidx *index,struct input query)
{
        FILE *f;
	int i,j;
	struct matches *temp;
        if(strlen(OUTPUT_FILE) == 0)
        {
                f = fopen("out.bed","w");
        }
        else
        {
                strcat(OUTPUT_FILE,".bed");
                f = fopen(OUTPUT_FILE,"w");
        }
/*      write each instance of M to a file*/
        if (f == NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	fprintf(f,"Database\tStart Pos\tEnd Pos\tQuery\tScore\n----------------------------------------------------------------\n");
        for(i = 0; i < qsc; i++)
        {
                for(j = 0; j < isc; j++)
                {
			temp = (struct matches *) malloc(sizeof(struct matches));
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

/*char *getSequenceAlignment(int i,int j,struct FMidx *index, struct input query,struct matches *match)
{
	int z,q,miss,nLcount,temp;
	char *sequence;
	if(match->traceLength > 60)
	{
		nLcount = (match->traceLength/60) + 3;
		sequence = (char *) malloc((nLcount+match->traceLength)*sizeof(char));
	}
	else
	{
		sequence = (char *) malloc(match->traceLength*sizeof(char));
	}
	temp = 0;
	for(q = 0; q < nLcount; q++)
	{
		miss = 0;
		for(z = 0; z < match->traceLength; z++)
		{
			if(z < 60)
			{
				if(match->tbi[z] == 1 || match->tbi[z] == 3)
				{
					sequence[z] = index[j].sequence[index[j].SA[match->low]+(z-1-miss)];
					printf("%c\n",index[j].sequence[index[j].SA[match->low]+(z-1-miss)]);
				}
				else if(match->tbi[z] == 2)
				{
					sequence[z] = '-';
					miss++;
				}
			}
			else
			{
				sequence[z] = '\n';
				temp = z+1+temp;
				z = z + match->traceLength;
			}
		}
		miss = 0;
		for(z = 0; z < match->traceLength; z++)
                {
			if(z < 60)
                        {
                                if(match->tbi[z] == 1)
                                {
					if(index[j].sequence[index[j].SA[match->low]+(z-1-miss)] == query.sequence[i][z-miss])
					{
						sequence[z+temp] = '|';
					}
					else
					{
						sequence[z+temp] = ' ';
					}
                                }
				else if(match->tbi[z] == 3 || match->tbi[z] == 2)
				{
					sequence[z+temp] = ' ';
				}
                        }
                        else
                        {
                                sequence[z] = '\n';
                                temp = z+1+temp;
                                z = z + match->traceLength;
                        }
		}
		miss = 0;
                for(z = 0; z < match->traceLength; z++)
                {
                        if(z < 60)
                        {
                                if(match->tbi[z] == 1 || match->tbi[z] == 2)
                                {
                                        sequence[z+temp] = query.sequence[i][z-miss];
                                }
                                else if(match->tbi[z] == 3)
                                {
                                        sequence[z] = '-';
                                        miss++;
                                }
                        }
                        else
                        {
                                sequence[z] = '\n';
				sequence[z+1] = '\n';
                                z = z + match->traceLength;
                        }
                }
	}
	return sequence;
}*/

/*				*
 *	PRINTING FUNCTIONS	*
 * 				*/

/*
 
	printInResults

	out:

	qsc:

	isc:

	index:

	query:

*/ 
void printInResults(struct results **out,int qsc,int isc,struct FMidx *index, struct input query)
{
	int i,j,z;
	int miss = 0;
	int missX = 0;
        int missY = 0;
	struct matches *temp;
	for(i = 0; i < qsc; i++)
	{
		for(j = 0; j < isc; j++)
		{
			temp = (struct matches *) malloc(sizeof(struct matches));
			temp = out[i][j].match; /*causing seg fault*/
			printf("\nindex: %s\nquery: %s\n\n-----------------------------------------------\n\n",index[j].desc,query.name[i]);
			if(temp!=NULL)
			{
				while(temp!=NULL)
				{
/*					printf("tb:%s\n",temp->tb);*/
					printf("Score: %d\n\nIndex\t",temp->score);
		/*			char *sequence = getSequenceAlignment(i,j,index,query,temp);
					printf("%s",sequence);*/
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
					missX = 0;
					missY = 0;
					for(z = 0; z < temp->traceLength; z++)
                                        {
						if(z % 59 == 0 && z != 0)
                                                {
                                                        printf("\n");
                                                }
                                                else
                                                {
							if(temp->tb[z] == 'X')
							{
								printf(" ");
								missX++;
							}
							else if(temp->tb[z] == 'Y')
                                                        {
                                                                printf(" ");
                                                                missY++;
                                                        }
							else if(temp->tb[z] == 'U')
							{
								printf(" ");
							}
                                                        else if(temp->tb[z] == 'M')
							{
                                                                printf("|");
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

					printf("\n\nStart: %d\nEnd: %d\n\n-----------------------------------------------\n\n\n\n",index[j].SA[temp->low],index[j].SA[temp->low]+temp->traceLength-1);
					temp = temp->next;
				}
			}
			else
			{
				printf("No Matches Found\n\n-----------------------------------------------\n\n\n\n");
			}
		}
	}	
}

/*

	printKLs

	prints the calculated K and L values from each input/index combination

	out: 2D array of structures that contain search output information

	qsc: integer variable containing the number of query sequences

	isc: integer variable containing the number of index sequences

*/

void printKLs(struct results **out, int qsc, int isc)
{
	int i,j;
	struct matches *temp;
	for(i = 0; i < qsc; i ++)
	{
		for(j = 0; j < isc; j++)
		{
			printf("index: %d\tquery: %d\n",j+1,i+1);
			temp = (struct matches *) malloc(sizeof(struct matches));
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

/*				*
 * 	MAIN FUNCTION		*
 * 				*/

int main(int argc, char *argv[])
{
	int QseqCount = 0;
	int IseqCount = 0;
	int ***D; 
        struct results **out;
	struct input query = manageInputs(argv,argc,&QseqCount);
        struct FMidx *index = getIndex();
	IseqCount = getCount();
	D = calculateD(index,IseqCount,query,QseqCount);
	/*Search*/
/*	struct output **out = search(query,QseqCount,index,IseqCount);*/
	out = inexactSearch(query,QseqCount,index,IseqCount,D);
	/*Output*/
	outputToFile(out,QseqCount,IseqCount,index,query);
	if(verbose)
	{
		printKLs(out,QseqCount,IseqCount);
	}
	if(!silent)
	{
		printInResults(out,QseqCount,IseqCount,index,query);
	}
	return 0;
}
