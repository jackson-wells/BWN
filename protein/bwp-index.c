#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "bwp-index.h"
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <getopt.h>

bool verbose = false;
int MAX_LINE_LENGTH = 10000000;
char OUTPUT_FILE[] = "index.bwp";

/* 				*
 * 	PRINTING FUNCTIONS 	*
 * 				*/


/*
 
	printSuffixArray

	outputs the suffix array(s) calculated from the input file

	m: 2 dimensional structure array containing the sorted suffix array(s)

	seqCount: integer variable containing the number of seqeunces presint
                in the input file

        seqLength: an array of integers containing the lengths of each sequence
                in the input file

	
*/
void printSuffixArray(struct suffix **m,int seqCount, int *seqLength)
{
	int i,j;
	for(i= 0; i < seqCount;i++)
	{
		printf("\nIndex Sequence: %d\n",i);
		for(j=0; j < seqLength[i]; j++)
		{
			printf("%d\t%d\t%s\n",j,m[i][j].pos-1,m[i][j].string);
		}
	}
	printf("\n");
}

/*
 
	printBwt

	outputs the transform(s) calculated from the input file

	transform:  2 dimensional character array containing the
                burrows-wheeler transform of each sequence from the input file

	seqCount: integer variable containing the number of seqeunces presint
                in the input file

*/
void printBwt(char **transform,int seqCount)
{
	int i;
	for(i = 0;i < seqCount; i++)
    	{
        	printf("bwt:\n%s\n",transform[i]);
    	}
	printf("\n");
}


/* 					*
 * 	MEMORY ALLOCATION FUNCTIONS 	*
 * 					*/


/*
 
	initializeSuffixArray

	allocates memory for the variables contained in the suffix
	structure object

	m:  structure array to contain the suffix array

	seqLength: an array of integers containing the lengths of each sequence
                in the input file


*/
void initializeSuffixArray(struct suffix *m,int seqLength)
{
	int j;
	for(j = 0; j < seqLength; j++)
        {
                m[j].string = (char *) malloc(seqLength * sizeof(char));
    	        m[j].pos = j+1;
                assert(m[j].string != 0);
        }
}

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
struct input initializeInputStruct( int seqCount, int *seqLength)
{
	struct input query;
    	int i;
    	int n = 100;
    	query.length = (int *) malloc(seqCount * sizeof(int));
    	query.name = (char **) malloc(seqCount * sizeof(char *));
	query.sequence = (char **) malloc(seqCount * sizeof(char *));
	query.reverse = (char **) malloc(seqCount * sizeof(char *));
	for(i = 0; i < seqCount;i++)
	{
		query.name[i] = (char *) malloc(n*(sizeof(char)));
	        query.sequence[i] = (char *) malloc(seqLength[i] * sizeof(char));
		query.reverse[i] = (char *) malloc(seqLength[i] * sizeof(char));
	}
    	return query;
}



/* 					*
 * 	MEMORY CLEARING FUNCTIONS 	*
 * 					*/


/*
 
	deleteInputStruct

	manages the removal of the input structure variable from memory

	query: structure variable containing file information

        seqCount: integer variable containing the number of seqeunces presint
                in the input file

*/
void deleteInputStruct(struct input query,int seqCount)
{
    	int i;
    	for(i = 0; i < seqCount; i++)
    	{
/*    	    free(query.name[i]);
    	    free(query.sequence[i]);*/
    	    query.length[i] = 0;
    	}
}


/*
 
	deleteSuffixArray

	manages the removal of the suffix structure variable from memory

	m: 2 dimensional structure array containing the sorted suffix array(s)

        seqCount: integer variable containing the number of seqeunces presint
                in the input file

        seqLength: an array of integers containing the lengths of each sequence
                in the input file


*/
void deleteSuffixArray(struct suffix **m,int seqCount,int *seqLength)
{
	int i;
	freeSuffixArray(m,seqCount,seqLength);
	for(i = 0; i< seqCount; i++)
	{
	/*	free(m[i]);*/
	}
}

/*
 
	freeSuffixArray

	manages the removal of the varaibles contained in the suffix
	structure object from memory

	m: 2 dimensional structure array containing the sorted suffix array(s)

        seqCount: integer variable containing the number of seqeunces presint
                in the input file

        seqLength: an array of integers containing the lengths of each sequence
                in the input file

*/
void freeSuffixArray(struct suffix **m,int seqCount,int *seqLength)
{
	int i,j;
	for(i = 0; i < seqCount; i++)
	{
		for(j = 0;j< seqLength[i];j++)
		{
			if(m[i][j].string != 0)
			{
			/*	free(m[i][j].string); */	/* free the space on the heap */
				m[i][j].string = 0;   	/* make it point to null */
			}
			m[i][j].pos = 0;
		}
	}
}



/* 				*
 * 	HELPER FUNCTIONS 	*
 * 				*/

/*
  
 	fileExists
	
	verifies the input file exists and is accessible

	temp: string variable containing the name of the input file

	returns: 
		0 if the file exists and is accessible
		
		1 if the file does not exist or is unaccessible

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

/*  
 
	seqCount

    	finds the number of sequences in the input file
	
	fileName: string variable containing the name of the input file

	returns: integer corresponding to the number of sequences contained
		in the input file

*/
int getSeqCount(char *fileName)	/*reutrns number of sequences present in a multi-fasta file*/
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
	/*free(temp);*/
	return lineCount;
}

/*  

	seqLength

    	gets the lengths of query sequences from the input file

	fileName: string variable containing the name of the input file	
	
	seqCount: integer variable containing the number of sequences 
		contained in the input file

	returns: an array of integers

*/
int *getLength(char *fileName,int seqCount) /*returns an array of sequence lengths*/
{
	int i = 0;
        FILE *file = fopen(fileName,"r");
        char *seq = (char *) malloc(MAX_LINE_LENGTH * sizeof(char));
        int *seqLength = (int *) malloc(seqCount * sizeof(int));;
        while(fgets(seq,MAX_LINE_LENGTH-1,file) != NULL)
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

/*  

 	read_fasta

   	stores fasta file information into memory, adds "$" character
	to sequence and reverse sequence
	
	fileName: string variable containing the name of the file input
		to be indexed

	query: structure variable to contain file information

*/
void read_fasta(char *fileName, struct input *query)
{
    	char *temp = (char *) malloc(MAX_LINE_LENGTH * sizeof(char));
	char *rev = (char *) malloc(MAX_LINE_LENGTH * sizeof(char)); 
    	int i = 0;
	FILE *file = fopen(fileName,"r");
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
            	else    /*if line contains an amino acid sequence*/
            	{
                    	strtok(temp,"\n"); /*strings read from file have extra \n added by file read*/
			rev = reverse(temp);
			strcat(rev,"$");
                    	strcat(temp,"$");
                    	strncpy(query->sequence[i],temp,strlen(temp));    /*saving string in memory*/
			strncpy(query->reverse[i],rev,strlen(rev));
                    	query->length[i] = strlen(temp);
                    	i++;
			memset(temp,0,strlen(temp));
			memset(rev,0,strlen(rev));
            	}
    	}
	fflush(file);
	fclose(file);
/*	free(temp);
	free(rev);*/
}

void formatFasta(char *fileName)
{
	int i;
	int lineCount = getLineCount(fileName);
	char *temp = (char *) malloc(MAX_LINE_LENGTH * sizeof(char));
	char **fileInfo = (char **) malloc(lineCount * sizeof(char *));
	FILE *file;
	for(i = 0; i < lineCount; i++)
	{
		fileInfo[i] = (char *) malloc(MAX_LINE_LENGTH * sizeof(char)); 
	}
	i = 0;
	file = fopen(fileName,"r");
	while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)   /*loops through each line of a file*/
        {
		temp[strcspn(temp, "\n")] = 0;
		temp[strcspn(temp, "\r")] = 0;
		strcpy(fileInfo[i],temp);
		i++;
	}
	fclose(file);
	fflush(file);
	file = fopen(fileName,"w+");
	for(i = 0; i < lineCount; i++)
	{
		if(fileInfo[i][0] == '>' && i != 0)
		{
			fprintf(file,"\n%s\n",fileInfo[i]);
		}
		else if(fileInfo[i][0] == '>' && i == 0)
		{
			fprintf(file,"%s\n",fileInfo[i]);
		}
		else
		{
			fprintf(file,"%s",fileInfo[i]);
		}
	}
	fflush(file);
	fclose(file);
	for(i = 0; i < lineCount; i++)
        {
		free(fileInfo[i]);
	}
	free(temp);
	free(fileInfo);
}

int getLineCount(char *fileName)
{
	char *temp = (char *) malloc(MAX_LINE_LENGTH * sizeof(char));
	int lineCount = 0;
        FILE *file = fopen(fileName,"r");
        while(fgets(temp,MAX_LINE_LENGTH,file) != NULL)   /*loops through each line of a file*/
        {
		lineCount++;
	}
	fclose(file);
	return lineCount;
}


/*  

	charToEnd

	shifts a string left, maintaing the first element 
	and moving it to the back of the string

	input: string to be shifted 
	
*/
void charToEnd(char* input)
{			
	const int len = strlen(input);
    	if(len > 1)
    	{
        	const char first = input[0];
        	memmove(input,input+1,len-1);
        	input[len -1] = first;
    	}
}

/*	

	reverse
 	
	reverses the order of characters in a string

	str: string to be reversed

	returns: a temporary variable containing a reversed verison
		of the input string

*/
char* reverse(char *str)
{
	char *temp = (char *) malloc(strlen(str)*sizeof(char));
	int i;
	int idx = 0;
	for(i = strlen(str)-1; i > -1; i=i-1)
	{
		temp[idx] = str[i];
		idx++;
	}
	return temp;
}

/*
 	
	extensionExists

	iterates over the input file name to verify if a file
	extension exists

	temp: string variable containing the name of the input file

	returns:
		0 if the input file's name contains a file extension
		
		1 if the input file's name does not conatain a file extension 

*/
int extensionExists(char *temp)
{
        int i = 0;
        for(i = 0; i < (int)(strlen(temp)-1); i++)
        {
                if(temp[i] == '.')
                {
                        return 0;
                }
        }
        return 1;
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
        if(temp == 'A'){return 0;}
        else if(temp == 'C'){return 1;}
        else if(temp == 'D'){return 2;}
        else if(temp == 'E'){return 3;}
        else if(temp == 'F'){return 4;}
        else if(temp == 'G'){return 5;}
        else if(temp == 'H'){return 6;}
        else if(temp == 'I'){return 7;}
        else if(temp == 'K'){return 8;}
        else if(temp == 'L'){return 9;}
        else if(temp == 'M'){return 10;}
        else if(temp == 'N'){return 11;}
        else if(temp == 'P'){return 12;}
        else if(temp == 'Q'){return 13;}
        else if(temp == 'R'){return 14;}
        else if(temp == 'S'){return 15;}
        else if(temp == 'T'){return 16;}
        else if(temp == 'V'){return 17;}
        else if(temp == 'W'){return 18;}
        else if(temp == 'Y'){return 19;}
        else
        {
                printf("Invalid character '%c'\n",temp);
                exit(1);
        }
}

/* 					*
 * 	INPUT HANDLING FUNCTIONS 	*
 * 					*/


/*
 
	manageInputs

	handles command-line input from user

	argc: integer variable containing the number of arguements passed in

	argv: 2d character array containing the arguement strings

	seqCount: integer array to hold the number of sequences to be processed
		in the run of the program

	returns: a structure variable containing input information


*/
struct input manageInputs(int argc, char *argv[],int *seqCount) 
{
	struct input query;
  	int c;
	if(argc <= 1)	/*no arguements supplied*/
	{
		printf("Please provide the necessary options\n\nuse -h for usage statement\n");
		exit(0);
	}
    	opterr = 0;
    	while ((c = getopt (argc, argv, "hvf:s:m:o:")) != -1) /*options must be added here to be recognized, options followed by : take in a parameter*/
    	{	
	        switch (c)
	        {
			printf("hello?\n");
	            	case 'h':
	                	printf("\nBurrows Wheeler Nucleotide Alligner\n\nUsage: \"bwp-index <options>\"\n\nOptions:\n\n-f\t\tFor input of a fasta fileas a query\n-s\t\tFor input of a string as a query\n-h\t\tFor this usage statement\n-m\t\tTo designate maximum sequence length according to character count\n-o\t\tTo designate the output file name\n-v\t\tTo produce verbose output\n\n");
	                	exit(0);
	
	            	case 'f':
				if(fileExists(optarg))
                                {
		                	handleF(&query,optarg);
		                	*seqCount = getSeqCount(optarg);
				}
				else
				{
                                        printf("%s is not a valid file, exiting\n",optarg);
					exit(0);
                                }
		                break;
	
		        case 's' :
		                handleS(&query,optarg);
		                *seqCount = 1;
		                break;
			case 'v' :
                                verbose = true;
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
	int seqCount = 1;
	int *seqLength = (int *) malloc(sizeof(int));
	char *rev = reverse(sequence);
    	printf("\nHere is the string you entered:\n%s\n\n",sequence);
    	seqLength[0] = strlen(sequence) + 1;
	strcat(rev,"$");
    	strcat(sequence,"$");
    	*query = initializeInputStruct(seqCount,seqLength);
    	query->length[0] = seqLength[0];
	strcpy(query->reverse[0],rev);
    	strcpy(query->sequence[0],sequence);
    	strcpy(query->name[0],"Command-line Input");
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
    	int *seqLength = getLength(fileName,seqCount);
    	*query = initializeInputStruct(seqCount,seqLength);
	formatFasta(fileName);
    	read_fasta(fileName, query);
}



/* 					*
 *  	SUFFIX ARRAY FUNCTIONS 		*
 *  					*/


/*
 
	newSuffixArray

	allocates memory for a structure array and the variables contained
	in said structure array

	seqLength: an integer variable containing the lenght of the current
		input sequence

	returns: an initialized and empty structure array

*/
struct suffix *newSuffixArray(int seqLength)
{
        struct suffix *m = (struct suffix *) malloc(seqLength * sizeof(struct suffix));
        assert(m != 0);
        initializeSuffixArray(m,seqLength);
        return m;
}

/*
 
	buildSuffixArray

	manages memory declaration and sorting of a SA
	
	sequence: a character array containing an input sequence

	returns: a structure array containing a sorted suffix array


*/
struct suffix *buildSuffixArray(char *sequence)
{
    	struct suffix *SA = populateSuffixArray(sequence);
	int length = strlen(sequence);
/*	printf("here4.17\n");
	printf("len:%d\n",strlen(sequence));*/
    	mergeSort(0,length-1,SA);
/*	printf("here4.19\n");*/
    	return SA;
}


/*
 
	populateSuffixArray

	puts a sequence into a structre array, rotated each time

	sequence: a character array containing an input sequence

        returns: a structure array containing an unsorted suffix array


*/
struct suffix *populateSuffixArray(char *sequence)
{
	int length = strlen(sequence);
	struct suffix *SA = newSuffixArray(length);
	int j;
        char *tempSeq = sequence;
	length = strlen(sequence);
        for(j = 0; j < length; j++)
        {
                SA[j].pos = j + 1;
               	if(j == 0)
               	{
               		strcpy(SA[j].string,tempSeq);
               	}
               	if(j > 0)
               	{
               	    	charToEnd(tempSeq);
               	    	strcpy(SA[j].string,tempSeq);
               	}
        }
	return SA;
}

/* 				*
 * 	SORTING FUNCTIONS 	*
 * 				*/


/*
 
	mergeSort
	
	manages the recursive sorting of a suffix array

	low: an integer value initially holding the value of 0
	
	high: an integer value initially holding the 1 minus the length
		of the current sequence

	m: structure array containing the unsorted suffix array


*/
void mergeSort(int low, int high,struct suffix *m)
{
	if(low<high)
 	{
        	int mid = low+(high - low)/2;
        	mergeSort(low,mid,m);
        	mergeSort(mid+1,high,m);
        	Merge(low,mid,high,m);
    	}
}


/*

	Merge

	merges two structure arrays containing sequence and position
	information

	low: an integer value initially holding the value of 0

	mid: an integer value holding the middle location between the low
		and high variables

	high: an integer value initially holding the 1 minus the length
                of the current sequence

	m: structure array containing the unsorted suffix array


 */
void Merge(int low,int mid, int high, struct suffix *m)
{
	int j = 0;
	int k = low;
    	int nL= mid-low+1;
    	int nR= high-mid;
	int seqLength = strlen(m[0].string);
	int i = 0;
	struct suffix *tempL;
	struct suffix *tempR;
	tempL = (struct suffix *) malloc(sizeof(struct suffix)*nL);
	tempR = (struct suffix *) malloc(sizeof(struct suffix)*nR);
	for(i = 0; i < nL; i++)
	{
		tempL[i].string = malloc(seqLength * sizeof(char)); /*Memory allocation for string in struct*/
		strcpy(tempL[i].string,m[low+i].string); 
		tempL[i].pos = m[low+i].pos;
	}
	for(i = 0; i < nR; i++)
        {
		tempR[i].string = malloc(seqLength * sizeof(char));
                strcpy(tempR[i].string,m[mid+i+1].string);
                tempR[i].pos = m[mid+i+1].pos;
        }
	i = 0;
	while (i < nL && j < nR)
    	{
        	if (strcmp(tempL[i].string,tempR[j].string) <= 0)
        	{
        	    	strcpy(m[k].string,tempL[i].string);
			m[k].pos = tempL[i].pos;
        	    	i++;
        	}
        	else
        	{
			strcpy(m[k].string,tempR[j].string);
                        m[k].pos = tempR[j].pos;
        	    	j++;
        	}
        	k++;
    	}
	while(i < nL)
    	{
		strcpy(m[k].string,tempL[i].string);
                m[k].pos = tempL[i].pos;
        	i++;
        	k++;
    	}
	while (j < nR)
    	{
		strcpy(m[k].string,tempR[j].string);
                m[k].pos = tempR[j].pos;
        	j++;
        	k++;
    	}
}

/* 						*
 * 	TRANSFORM CALCULATION FUNCTIONS 	*
 * 						*/


/*
 
 	bwt

	calculates a transform sequence from a suffix array

	m: 2 dimensional structure array containing the sorted suffix array

	seqCount: integer variable containing the number of seqeunces presint
                in the input file

	seqLength: an array of integers containing the lengths of each sequence
		in the input file
	
	returns: a character array containing the calculated burrows-wheeler
		transform
  
 */
char **bwt(struct suffix **m, int seqCount,int *seqLength)
{
	int i,x;
	char **temp = (char **) malloc(seqCount * sizeof(char *));
	for(x = 0; x < seqCount; x++)
	{
		temp[x] = (char *) malloc(seqLength[x] * sizeof(char));
		for(i =0; i < seqLength[x]; i++)
		{
			/*printf("%c",m[x][i].string[seqLength[x]-1]);*/
			temp[x][i] = m[x][i].string[seqLength[x]-1];	/*gets last element of char* in each structure element*/
		}
	}
	return temp;
}



/* 						*
 * 	INTERVAL CALCULATION FUNCTIONS 		*
 * 						*/


/*
	calculateInterval

	manages the calculation of the O and C values

	transform: 2 dimensional character array containing the
                burrows-wheeler transform of each sequence from the input file
	
	seqCount: integer variable containing the number of seqeunces presint
                in the input file

	revTransform: 2 dimensional character array containing the
                reverse burrows-wheeler transform of each sequence from
                the input file

	returns: a structure array containing the calculated O and C values for
		each sequence in the input file

*/
struct FMidx *calculateInterval(char **transform, int seqCount, char **revTransform)
{
    	struct FMidx *index = (struct FMidx *) malloc(seqCount * sizeof(struct FMidx));
    	int i,j,length;
    	for(i = 0; i < seqCount; i++)
    	{
		index[i].O = (int **) malloc(20*sizeof(int **));
		index[i].R = (int **) malloc(20*sizeof(int **));
        	for(j = 0; j < 20; j++) /*looping thru each char*/
        	{
			length = strlen(transform[i]);
			index[i].O[j] = (int *) malloc(length*sizeof(int));
			index[i].R[j] = (int *) malloc(length*sizeof(int));
			index[i].R[j] = calculateO(revTransform[i],j);
           		index[i].O[j] = calculateO(transform[i],j);
			if(j == 0)				
/*Base Cases covered to decrease runtime*/
			{
				index[i].C[j] = 0;
			}
			else
			{
            			index[i].C[j] = calculateC(transform[i],j);
			}
        	}
    	}
    	return index;
} 


/*
 
	calculateO
	
	computes O values for a single character and sequence

	sequence: a string variable containing the transform of a sequence
                from the input file

	letterValue: an integer variable containing the number of the character
                being counted

	returns: an array of integers increasing from zero as the character
		is found in the transform, excluding "$"

			
*/
int *calculateO(char *sequence,int letterValue)
{
	int i;
        int count = 0;
	int length = strlen(sequence);
	int *value = (int *) malloc(length * sizeof(int)); 
        for(i = 0; i < length; i++)
        {
		if(sequence[i] != '$')
                {
                	if(baseMap(sequence[i]) == letterValue)
                	{
                	        count++;
                	}
		}
		value[i] = count;
        }
        return value;
}

/*
 
	calculateC

	computes C value for a single character
	
	sequence: a string variable containing the transform of a sequence
		from the input file

	letterValue: an integer variable containing the number of the character
		being counted
	
	returns: the number of occurrences of the character in the transform,
                excluding "$"

 */
int calculateC(char *sequence, int letterValue)
{
    	int i;
	int count = 0;
	int temp = strlen(sequence);
        for(i = 0; i < temp; i++)
        {
		if(sequence[i] != '$')
        	{	
			if(baseMap(sequence[i]) < letterValue)
           		{	
                		count++;
            		}
		}
        }
    	return count;
}

/* 				*
 * 	OUTPUT FUNCTIONS 	*
 * 				*/

/*
 
	intervalToFile
	
	writes the calculated sequence information to the specified file name

	index: structure array containing the O, R and C values calculated
		from the input file
	
	query: structure variable containing file information
	
	m: 2 dimensional structure array containing the sorted suffix array 
	
	transform: 2 dimensional character array containing the 
		burrows-wheeler transform of each sequence from the input file
	
	revTransform: 2 dimensional character array containing the
        	reverse burrows-wheeler transform of each sequence from
		the input file

	seqCount: integer variable containing the number of seqeunces presint
        	in the input file


*/
void intervalToFile(struct FMidx *index, int seqCount,struct suffix **m, char **transform, struct input query,char **revTransform)
{
	FILE *f;
	int i,j,z,q;
	if(strlen(OUTPUT_FILE) == 0) 
	{
		f = fopen("index.bwp","w");	/*change when moved to new infrastructure*/
	}
	else
	{
		if(extensionExists(OUTPUT_FILE))
		{
			strcat(OUTPUT_FILE,".bwp");
		}
		f = fopen(OUTPUT_FILE,"w");
	}
	/*write each instance of M to a file*/
	if (f == NULL)
	{
            	printf("Error opening file!\n");
        	exit(1);
	}
	fprintf(f,"n:%d\n\n",seqCount);
	for(i = 0; i<seqCount;i++)
        {
/*		query.sequence[i] = query.sequence[i]+2;*/
		fprintf(f,"d:%s\nl:%d\nq:%s\ns:",query.name[i],query.length[i],query.sequence[i]+1);
		for(q = 0; q < query.length[i]; q++)
		{
			fprintf(f,"%d ",m[i][q].pos);
		}
                fprintf(f,"\nt:%s\nf:%s\nc:%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",transform[i],revTransform[i],index[i].C[0],index[i].C[1],index[i].C[2],index[i].C[3],index[i].C[4],index[i].C[5],index[i].C[6],index[i].C[7],index[i].C[8],index[i].C[9],index[i].C[10],index[i].C[11],index[i].C[12],index[i].C[13],index[i].C[14],index[i].C[15],index[i].C[16],index[i].C[17],index[i].C[18],index[i].C[19]);
                for(j= 0; j < 20; j++)
                {
			fprintf(f,"o:");
                        for(z= 0; z < query.length[i]; z++)
                        {
                                fprintf(f,"%d ",index[i].O[j][z]);
                        }
                        fprintf(f,"\n");
                }
		for(j= 0; j < 20; j++)
                {
                        fprintf(f,"r:");
                        for(z= 0; z < query.length[i]; z++)
                        {
                                fprintf(f,"%d ",index[i].R[j][z]);
                        }
                        fprintf(f,"\n");
                }
                fprintf(f,"\n");
        }
	fclose(f);
}



/* 			*
 *  	MAIN FUNCTION 	*
 *  			*/


int main(int argc, char *argv[])
{
	int seqCount = 0;	/*will contain # of sequences, is written by manageInputs()*/
	char **revTransform;
	char **transform;
	struct FMidx *index;
	struct input query = manageInputs(argc,argv,&seqCount);
	struct suffix **m = (struct suffix **) malloc(seqCount*sizeof(struct suffix *));
	struct suffix **Rm = (struct suffix **) malloc(seqCount*sizeof(struct suffix *));
	int i;
	for(i = 0; i < seqCount; i++)
	{
		m[i] = buildSuffixArray(query.sequence[i]);
/*		printf("here4.1\n");*/
		Rm[i] = buildSuffixArray(query.reverse[i]);
/*		printf("here4.2\n");*/
	}
/*	printf("here5\n");*/
	if(verbose)
	{
		printSuffixArray(m,seqCount,query.length);
	}
	transform = bwt(m,seqCount,query.length);	/*contains transforms of all query fastas*/
/*	printf("here5.5\n");*/
	revTransform = bwt(Rm,seqCount,query.length);
/*	printf("here6\n");*/
	if(verbose)
	{
		printBwt(transform,seqCount);
	}
/*	printf("here7\n");*/
	index = calculateInterval(transform,seqCount,revTransform);
/*	printf("here8\n");*/
	intervalToFile(index,seqCount,m,transform,query,revTransform);
/*	printf("here9\n");
    	deleteSuffixArray(m,seqCount,query.length);
    	deleteInputStruct(query,seqCount);*/
    	return 0;
}
