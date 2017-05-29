#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "bwn-index.h"
#include <unistd.h>
#include <math.h>

int MAX_LINE_LENGTH = 10000000;
char OUTPUT_FILE[100];

/* PRINT FUNCTIONS */

void printSuffixArray(struct suffix **m,int seqCount, int *seqLength)
{
	int i,j;
	for(i= 0; i < seqCount;i++)
	{
		/*printf("\n%s\n",m[i][0].name);*/
		for(j=0; j < seqLength[i]; j++)
		{
			printf("%d\t%d\t%s\n",j+1,m[i][j].pos,m[i][j].string);
		}
	}
	printf("\n");
}

void printBwt(struct suffix **m,char **transform,int seqCount)
{
	int i;
	for(i = 0;i < seqCount; i++)
    	{
            	/*printf("sequence:\tbwt:\n%s\t\t%s\n\n",m[i][0].name,transform[i]);*/
        	printf("bwt:\n%s\n",transform[i]);
    	}
	printf("\n");
}

void printInt(struct FMidx *index,int seqCount,int *seqLength)
{
        int i,j,z;
        for(i = 0; i<seqCount;i++)
        {
		printf("C(\"$\") of seq %d\t%d\n",i+1,index[i].C[0]);
            	printf("C(\"A\") of seq %d\t%d\n",i+1,index[i].C[1]);
            	printf("C(\"C\") of seq %d\t%d\n",i+1,index[i].C[2]);
            	printf("C(\"G\") of seq %d\t%d\n",i+1,index[i].C[3]);
            	printf("C(\"T\") of seq %d\t%d\n\n",i+1,index[i].C[4]);
		printf("O():\n");
		for(j= 0; j < 5; j++)
		{
			for(z= 0; z < seqLength[i]; z++)
			{
            			printf("%d ",index[i].O[j][z]);
/*            			printf("0(\"C\") of seq %d\t%d\n",i+1,index[i].O[1]);
			        printf("0(\"G\") of seq %d\t%d\n",i+1,index[i].O[2]);
			        printf("0(\"T\") of seq %d\t%d\n\n",i+1,index[i].O[3]);*/
			}
			printf("\n");
		}
        }
}



/* MEMORY ALLOCATION FUNCTIONS */

struct suffix **newSuffixArray(int seqCount, int *seqLength)
{
	struct suffix **m = (struct suffix **) malloc(seqCount * sizeof(struct suffix *));
	int i;
	for(i = 0; i < seqCount; i++) /*initiallizing x arrays of structs*/
	{
		m[i] = (struct suffix *) malloc(seqLength[i] * sizeof(struct suffix));
        	assert(m != 0);
	}
        initializeSuffixArray(m,seqLength,seqCount); /*initializes each element of x arrays*/
	return m;
}

void initializeSuffixArray(struct suffix **m,int *seqLength,int seqCount)
{
	int i,j;
	for(i = 0; i < seqCount;i++)
	{
		for(j = 0; j < seqLength[i]; j++)
        	{
        	        m[i][j].string = (char *) malloc(seqLength[i] * sizeof(char));
    	                m[i][j].pos = j+1;
        	        assert(m[i][j].string != 0);
        	}
	}
}

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



/* MEMORY CLEARING FUNCTIONS */

void deleteInputStruct(struct input query,int seqCount)
{
    	int i;
    	for(i = 0; i < seqCount; i++)
    	{
    	    free(query.name[i]);
    	    free(query.sequence[i]);
    	    query.length[i] = 0;
    	}
}

void deleteSuffixArray(struct suffix **m,int seqCount,int *seqLength)
{
	freeSuffixArray(m,seqCount,seqLength);
	int i;
	for(i = 0; i< seqCount; i++)
	{
		free(m[i]);
	}
}

void freeSuffixArray(struct suffix **m,int seqCount,int *seqLength)
{
	int i,j;
	for(i = 0; i < seqCount; i++)
	{
		for(j = 0;j< seqLength[i];j++)
		{
			if(m[i][j].string != 0)
			{
				free(m[i][j].string); 	/* free the space on the heap */
				m[i][j].string = 0;   	/* make it point to null */
			}
			m[i][j].pos = 0;
		}
	}
}



/* HELPER FUNCTIONS */

/*  seqCount

    Purpose: Find the number of sequences in the query file
	Returns: integer
*/
int getSeqCount(char *fileName)	/*reutrns number of sequences present in a multi-fasta file*/
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

/*  seqLength

    Purpose: get the lengths of query sequences from file
	Returns: array of integers
*/
int *getLength(char *fileName,int seqCount,int charCount) /*returns an array of sequence lengths*/
{
        FILE *file = fopen(fileName,"r");
        char *seq = (char *) malloc(charCount * sizeof(char));
        int *seqLength = (int *) malloc(seqCount * sizeof(int));;
        int i = 0;
        while(fgets(seq,charCount-1,file) != NULL)
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

    Purpose: Store FASTA file information into memory
	Returns: Nothing, Variables passed by reference
*/
void read_fasta(char *fileName, struct input *query)
{
    	int charCount = MAX_LINE_LENGTH;
    	char *temp = (char *) malloc(charCount * sizeof(char));
	char *rev = (char *) malloc(charCount * sizeof(char));
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
			rev = reverse(temp,strlen(temp));
			strcat(rev,"$");
                    	strcat(temp,"$");
                    	strcpy(query->sequence[i],temp);    /*saving string in memory*/
			strcpy(query->reverse[i],rev);
                    	query->length[i] = strlen(temp);
                    	i++;
            	}
    	}
	fflush(file);
    	fclose(file);
	//free(temp);
	//free(rev);
}

/*  charToEnd

    Purpose: Move first element of string to the last element.
	Returns: "Rotated" string
*/
void charToEnd(char* input)	/*Takes in char* and puts the first character element*/
{				/*at the back of the array */
	const int len = strlen(input);
    	if(len > 1)
    	{
        	const char first = input[0];
        	memmove(input,input+1,len-1);
        	input[len -1] = first;
    	}
}

char* reverse(char *str,int length)
{
	char *temp = (char *) malloc(length*sizeof(char));
	int i;
	int idx = 0;
	for(i = length-1; i > -1; i=i-1)
	{
		temp[idx] = str[i];
		idx++;
	}
	return temp;
}

/* INPUT HANDLING FUNCTIONS */

struct input manageInputs(int argc, char *argv[],int *seqCount) /*handles initial query into the program*/
{
	struct input query;
  	int c;
	if(argc <= 1)	/*no arguements supplied*/
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
	                	printf("\nBurrows Wheeler Nucleotide Alligner\n\nUsage: \"bwn-index <options>\"\n\nOptions:\n\n-f\t\tFor query of a fasta file\n-s\t\tFor query of a string\n-h\t\tFor this usage statement\n-m\t\tTo designate maximum sequence length according to character count\n-o\t\tTo designate the output file name\n\n");
	                	exit(0);
	
	            	case 'f':
		                handleF(&query,optarg);
		                *seqCount = getSeqCount(optarg);
		                break;
	
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
			case 'o' :
				strcpy(OUTPUT_FILE,optarg);
				
        	}
    	}
    	return query;
}

void handleS(struct input *query, char *sequence)
{
    	printf("\nHere is the string you entered:\n%s\n\n",sequence);
    	int *seqLength = (int *) malloc(sizeof(int));
    	seqLength[0] = strlen(sequence) + 1;
    	int seqCount = 1;
	char *rev = reverse(sequence,strlen(sequence));
	strcat(rev,"$");
    	strcat(sequence,"$");
    	*query = initializeInputStruct(seqCount,seqLength);
    	query->length[0] = seqLength[0];
	strcpy(query->reverse[0],rev);
    	strcpy(query->sequence[0],sequence);
    	strcpy(query->name[0],"Command Line String Input");
}

void handleF(struct input *query,char *fileName)
{
    	int seqCount = getSeqCount(fileName);
    	int charCount = MAX_LINE_LENGTH;
    	int *seqLength = getLength(fileName,seqCount,charCount);
    	*query = initializeInputStruct(seqCount,seqLength);
    	read_fasta(fileName, query);
}



/* SUFFIX ARRAY FUNCTIONS */

struct suffix **buildSuffixArray(struct input query,int seqCount)
{
    	/*intialize array*/
    	struct suffix **m = newSuffixArray(seqCount,query.length);
    	/*place sequences in respective areas of memory*/
    	populateSuffixArray(m,query,seqCount);
    	/*sort*/
    	int i;
    	for(i = 0; i < seqCount; i++)
    	{
    		mergeSort(0,query.length[i]-1,m,i);
    	}
    	return m;
}

struct suffix **buildReverseSuffixArray(struct input query,int seqCount)
{
        /*intialize array*/
        struct suffix **Rm = newSuffixArray(seqCount,query.length);
        /*place sequences in respective areas of memory*/
        populateReverseSuffixArray(Rm,query,seqCount);
        /*sort*/
        int i;
        for(i = 0; i < seqCount; i++)
        {
                mergeSort(0,query.length[i]-1,Rm,i);
        }
        return Rm;
}

void populateReverseSuffixArray(struct suffix **m, struct input query, int seqCount)
{
        int i,j;
        for(i = 0; i < seqCount; i++)
        {
                char *tempSeq = query.reverse[i];
                for(j = 0; j < query.length[i]; j++)
                {
                        m[i][j].pos = j + 1;
                        if(j == 0)
                        {
                                strcpy(m[i][j].string,tempSeq);
                        }
                        if(j > 0)
                        {
                                charToEnd(tempSeq);
                                strcpy(m[i][j].string,tempSeq);
                        }
                }
        }
}

void populateSuffixArray(struct suffix **m, struct input query, int seqCount)
{
	int i,j;
    	for(i = 0; i < seqCount; i++)
    	{
        	char *tempSeq = query.sequence[i];
        	for(j = 0; j < query.length[i]; j++)
        	{
                	m[i][j].pos = j + 1;
                	if(j == 0)
                	{
                		strcpy(m[i][j].string,tempSeq);
                	}
                	if(j > 0)
                	{
                	    	charToEnd(tempSeq);
                	    	strcpy(m[i][j].string,tempSeq);
                	}
        	}
    	}
}

void mergeSort(int low, int high,struct suffix **m,int seqNumber)
{
	if(low<high)
 	{
        	int mid = low+(high - low)/2;
        	mergeSort(low,mid,m,seqNumber);
        	mergeSort(mid+1,high,m,seqNumber);
        	Merge(low,mid,high,m,seqNumber);
    	}
}

void Merge(int low,int mid, int high, struct suffix **m, int seqNumber)
{
    	int nL= mid-low+1;
    	int nR= high-mid;
	int seqLength = strlen(m[seqNumber][0].string);
	int i = 0;
	struct suffix tempL[nL];
	struct suffix tempR[nR];
	for(i = 0; i < nL; i++)
	{
		tempL[i].string = malloc(seqLength-1 * sizeof(char)); //Memory allocation for string in struct
		strcpy(tempL[i].string,m[seqNumber][low+i].string); 
		tempL[i].pos = m[seqNumber][low+i].pos;
	}
	for(i = 0; i < nR; i++)
        {
		tempR[i].string = malloc(seqLength-1 * sizeof(char));
                strcpy(tempR[i].string,m[seqNumber][mid+i+1].string);
                tempR[i].pos = m[seqNumber][mid+i+1].pos;
        }
	int j = 0;
	i = 0;
	int k = low;
	while (i < nL && j < nR)
    	{
        	if (strcmp(tempL[i].string,tempR[j].string) <= 0)
        	{
        	    	strcpy(m[seqNumber][k].string,tempL[i].string);
			m[seqNumber][k].pos = tempL[i].pos;
        	    	i++;
        	}
        	else
        	{
			strcpy(m[seqNumber][k].string,tempR[j].string);
                        m[seqNumber][k].pos = tempR[j].pos;
        	    	j++;
        	}
        	k++;
    	}
	while(i < nL)
    	{
		strcpy(m[seqNumber][k].string,tempL[i].string);
                m[seqNumber][k].pos = tempL[i].pos;
        	i++;
        	k++;
    	}
	while (j < nR)
    	{
		strcpy(m[seqNumber][k].string,tempR[j].string);
                m[seqNumber][k].pos = tempR[j].pos;
        	j++;
        	k++;
    	}
}

/* TRANSFORM CALCULATION FUNCTIONS */

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



/* INTERVAL CALCULATION FUNCTIONS */

struct FMidx *calculateInterval(char **transform, int *seqLength,int seqCount, char **revTransform)
{
    	struct FMidx *index = (struct FMidx *) malloc(seqCount * sizeof(struct FMidx));
    	int i,j,z;
    	for(i = 0; i < seqCount; i++)
    	{
		index[i].O = (int **) malloc(5*sizeof(int **));
		index[i].R = (int **) malloc(5*sizeof(int **));
        	for(j = 0; j < 5; j++) //looping thru each char ($,A,C,G,T)
        	{
			index[i].O[j] = (int *) malloc(seqLength[i]*sizeof(int));
			index[i].R[j] = (int *) malloc(seqLength[i]*sizeof(int));
			index[i].R[j] = calculateO(revTransform[i],seqLength[i],j);
           		index[i].O[j] = calculateO(transform[i],seqLength[i],j);
			index[i].C[5] = seqLength[i]-1;
//			index[i].C[j] = index[i].O[j][seqLength[i]-1];
			if(j == 0)
			{
				index[i].C[j] = 0;
			}
			else if(j == 1)				//Base Cases covered to decrease runtime
			{
				index[i].C[j] = 1;
			}
			else
			{
            			index[i].C[j] = calculateC(transform[i],seqLength[i],j);
			}
        	}
    	}
    	return index;
} 

int baseMap(char temp)
{
	if(temp == '$') return 0;
	else if(temp == 'A') return 1;
	else if(temp == 'C') return 2;
	else if(temp == 'G') return 3;
	else if(temp == 'T') return 4;
}

int *calculateO(char *sequence,int seqLength,int letterValue)
{
	int i;
        int count = 0;
	int *temp = (int *) malloc(seqLength * sizeof(int)); 
        for(i = 0; i < seqLength; i++)
        {
                if(baseMap(sequence[i]) == letterValue)
                {
                        count++;
                }
		temp[i] = count;
        }
        return temp;
}

int calculateC(char *sequence, int seqLength,int letterValue)
{
    	int i;
	int count = 0;
        for(i = 0; i < seqLength; i++)
        {
        	if(baseMap(sequence[i]) < letterValue)
           	{
                	count++;
            	}
        }
    	return count;
}



/* OUTPUT FUNCTIONS */

void intervalToFile(struct FMidx *index, int seqCount,struct suffix **m, char **transform, struct input query)
{
	FILE *f;
	if(strlen(OUTPUT_FILE) == 0) 
	{
		f = fopen("/nfs0/Hendrix_Lab/bwp/nucl/index.bwn","w");	//change when moved to new infrastructure
	}
	else
	{
		strcat(OUTPUT_FILE,".bwn");
		f = fopen("/nfs0/Hendrix_Lab/bwp/nucl/OUTPUT_FILE","w");
	}
	/*write each instance of M to a file*/
	if (f == NULL)
	{
            	printf("Error opening file!\n");
        	exit(1);
	}
	int i,j,z,q;
	fprintf(f,"n:%d\n\n",seqCount);
	for(i = 0; i<seqCount;i++)
        {
		fprintf(f,"d:%s\nl:%d\ns:",query.name[i],query.length[i]);
		for(q = 0; q < query.length[i]; q++)
		{
			fprintf(f,"%d ",m[i][q].pos);
		}
                fprintf(f,"\nt:%s\nc:%d %d %d %d %d %d\n",transform[i],index[i].C[0],index[i].C[1],index[i].C[2],index[i].C[3],index[i].C[4],index[i].C[5]);
                for(j= 0; j < 5; j++)
                {
			fprintf(f,"o:");
                        for(z= 0; z < query.length[i]; z++)
                        {
                                fprintf(f,"%d ",index[i].O[j][z]);
                        }
                        fprintf(f,"\n");
                }
		for(j= 0; j < 5; j++)
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



/* MAIN FUNCTION */

int main(int argc, char *argv[])
{
	int seqCount = 0;	/*will contain # of sequences, is written by manageInputs()*/
	struct input query = manageInputs(argc,argv,&seqCount);
	struct suffix **m = buildSuffixArray(query,seqCount);
	struct suffix **Rm = buildReverseSuffixArray(query,seqCount);
	//printSuffixArray(m,seqCount,query.length);
	char **transform = bwt(m,seqCount,query.length);	/*contains transforms of all query fastas*/
	char **revTransform = bwt(Rm,seqCount,query.length);
	//printBwt(m,transform,seqCount);
	struct FMidx *index = calculateInterval(transform,query.length,seqCount,revTransform);
	//printInt(index,seqCount,query.length);
	intervalToFile(index,seqCount,m,transform,query);
    	deleteSuffixArray(m,seqCount,query.length);
    	deleteInputStruct(query,seqCount);
    	return 0;
}
