#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "batN-index.h"
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

void printInt(struct index *FMidx,int seqCount,int *seqLength)
{
        int i,j,z;
        for(i = 0; i<seqCount;i++)
        {
		printf("C(\"$\") of seq %d\t%d\n",i+1,FMidx[i].C[0]);
            	printf("C(\"A\") of seq %d\t%d\n",i+1,FMidx[i].C[1]);
            	printf("C(\"C\") of seq %d\t%d\n",i+1,FMidx[i].C[2]);
            	printf("C(\"G\") of seq %d\t%d\n",i+1,FMidx[i].C[3]);
            	printf("C(\"T\") of seq %d\t%d\n\n",i+1,FMidx[i].C[4]);
		printf("O():\n");
		for(j= 0; j < 5; j++)
		{
			for(z= 0; z < seqLength[i]; z++)
			{
            			printf("%d ",FMidx[i].O[j][z]);
/*            			printf("0(\"C\") of seq %d\t%d\n",i+1,FMidx[i].O[1]);
			        printf("0(\"G\") of seq %d\t%d\n",i+1,FMidx[i].O[2]);
			        printf("0(\"T\") of seq %d\t%d\n\n",i+1,FMidx[i].O[3]);*/
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

struct input_data initializeInputStruct( int seqCount, int *seqLength)
{
	struct input_data input;
    	int i;
    	int n = 100;
    	input.length = (int *) malloc(seqCount * sizeof(int));
    	input.name = (char **) malloc(seqCount * sizeof(char *));
	input.sequence = (char **) malloc(seqCount * sizeof(char *));
	input.reverse = (char **) malloc(seqCount * sizeof(char *));
	for(i = 0; i < seqCount;i++)
	{
		input.name[i] = (char *) malloc(n*(sizeof(char)));
	        input.sequence[i] = (char *) malloc(seqLength[i] * sizeof(char));
		input.reverse[i] = (char *) malloc(seqLength[i] * sizeof(char));
	}
    	return input;
}



/* MEMORY CLEARING FUNCTIONS */

void deleteInputStruct(struct input_data input,int seqCount)
{
    	int i;
    	for(i = 0; i < seqCount; i++)
    	{
    	    free(input.name[i]);
    	    free(input.sequence[i]);
    	    input.length[i] = 0;
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

    Purpose: Find the number of sequences in the input file
	Returns: integer
*/
int seqCount(char *fileName)	/*reutrns number of sequences present in a multi-fasta file*/
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

    Purpose: get the lengths of input sequences from file
	Returns: array of integers
*/
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

    Purpose: Store FASTA file information into memory
	Returns: Nothing, Variables passed by reference
*/
void read_fasta(char *fileName, struct input_data *input)
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
                    	strcpy(input->name[i],temp);
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
                    	strcpy(input->sequence[i],temp);    /*saving string in memory*/
			strcpy(input->reverse[i],rev);
                    	input->length[i] = strlen(temp);
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

struct input_data manageInputs(int argc, char *argv[],int *sCount) /*handles initial input into the program*/
{
	struct input_data input;
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
	                	printf("\nBurrows Wheeler Protein Alligner\n\nUsage: \"batN-index <options>\"\n\nOptions:\n\n-f\t\tFor input of a fasta file\n-s\t\tFor input of a string\n-h\t\tFor this usage statement\n-m\t\tTo designate maximum sequence length according to character count\n-o\t\tTo designate the output file name\n\n");
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

void handleS(struct input_data *input, char *sequence)
{
    	printf("\nHere is the string you entered:\n%s\n\n",sequence);
    	int *seqLength = (int *) malloc(sizeof(int));
    	seqLength[0] = strlen(sequence) + 1;
    	int seqC = 1;
	char *rev = reverse(sequence,strlen(sequence));
	strcat(rev,"$");
    	strcat(sequence,"$");
    	*input = initializeInputStruct(seqC,seqLength);
    	input->length[0] = seqLength[0];
	strcpy(input->reverse[0],rev);
    	strcpy(input->sequence[0],sequence);
    	strcpy(input->name[0],"Command Line String Input");
}

void handleF(struct input_data *input,char *fileName)
{
    	int sCount = seqCount(fileName);
    	int cCount = MAX_LINE_LENGTH;
    	int *sLength = seqLength(fileName,sCount,cCount);
    	*input = initializeInputStruct(sCount,sLength);
    	read_fasta(fileName, input);
}



/* SUFFIX ARRAY FUNCTIONS */

struct suffix **buildSuffixArray(struct input_data input,int seqCount)
{
    	/*intialize array*/
    	struct suffix **m = newSuffixArray(seqCount,input.length);
    	/*place sequences in respective areas of memory*/
    	populateSuffixArray(m,input,seqCount);
    	/*sort*/
    	int i;
    	for(i = 0; i < seqCount; i++)
    	{
    		mergeSort(0,input.length[i]-1,m,i);
    	}
    	return m;
}

struct suffix **buildReverseSuffixArray(struct input_data input,int seqCount)
{
        /*intialize array*/
        struct suffix **Rm = newSuffixArray(seqCount,input.length);
        /*place sequences in respective areas of memory*/
        populateReverseSuffixArray(Rm,input,seqCount);
        /*sort*/
        int i;
        for(i = 0; i < seqCount; i++)
        {
                mergeSort(0,input.length[i]-1,Rm,i);
        }
        return Rm;
}

void populateReverseSuffixArray(struct suffix **m, struct input_data input, int seqCount)
{
        int i,j;
        for(i = 0; i < seqCount; i++)
        {
                char *tempSeq = input.reverse[i];
                for(j = 0; j < input.length[i]; j++)
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

void populateSuffixArray(struct suffix **m, struct input_data input, int seqCount)
{
	int i,j;
    	for(i = 0; i < seqCount; i++)
    	{
        	char *tempSeq = input.sequence[i];
        	for(j = 0; j < input.length[i]; j++)
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

struct index *calculateInterval(char **transform, int *seqLength,int seqCount, char **revTransform)
{
    	struct index *FMidx = (struct index *) malloc(seqCount * sizeof(struct index));
    	int i,j,z;
    	for(i = 0; i < seqCount; i++)
    	{
		FMidx[i].O = (int **) malloc(5*sizeof(int **));
		FMidx[i].R = (int **) malloc(5*sizeof(int **));
        	for(j = 0; j < 5; j++) //looping thru each char ($,A,C,G,T)
        	{
			FMidx[i].O[j] = (int *) malloc(seqLength[i]*sizeof(int));
			FMidx[i].R[j] = (int *) malloc(seqLength[i]*sizeof(int));
			FMidx[i].R[j] = calculateO(revTransform[i],seqLength[i],j);
           		FMidx[i].O[j] = calculateO(transform[i],seqLength[i],j);
			FMidx[i].C[5] = seqLength[i]-1;
//			FMidx[i].C[j] = FMidx[i].O[j][seqLength[i]-1];
			if(j == 0)
			{
				FMidx[i].C[j] = 0;
			}
			else if(j == 1)				//Base Cases covered to decrease runtime
			{
				FMidx[i].C[j] = 1;
			}
			else
			{
            			FMidx[i].C[j] = calculateC(transform[i],seqLength[i],j);
			}
        	}
    	}
    	return FMidx;
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

void intervalToFile(struct index *FMidx, int seqCount,struct suffix **m, char **transform, struct input_data input)
{
	FILE *f;
	if(strlen(OUTPUT_FILE) == 0) 
	{
		f = fopen("/nfs0/Hendrix_Lab/bwp/nucl/index.batN","w");	
	}
	else
	{
		strcat(OUTPUT_FILE,".batN");
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
		fprintf(f,"d:%s\nl:%d\ns:",input.name[i],input.length[i]);
		for(q = 0; q < input.length[i]; q++)
		{
			fprintf(f,"%d ",m[i][q].pos);
		}
                fprintf(f,"\nt:%s\nc:%d %d %d %d %d %d\n",transform[i],FMidx[i].C[0],FMidx[i].C[1],FMidx[i].C[2],FMidx[i].C[3],FMidx[i].C[4],FMidx[i].C[5]);
                for(j= 0; j < 5; j++)
                {
			fprintf(f,"o:");
                        for(z= 0; z < input.length[i]; z++)
                        {
                                fprintf(f,"%d ",FMidx[i].O[j][z]);
                        }
                        fprintf(f,"\n");
                }
		for(j= 0; j < 5; j++)
                {
                        fprintf(f,"r:");
                        for(z= 0; z < input.length[i]; z++)
                        {
                                fprintf(f,"%d ",FMidx[i].R[j][z]);
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
	struct input_data input = manageInputs(argc,argv,&seqCount);
	struct suffix **m = buildSuffixArray(input,seqCount);
	struct suffix **Rm = buildReverseSuffixArray(input,seqCount);
	//printSuffixArray(m,seqCount,input.length);
	char **transform = bwt(m,seqCount,input.length);	/*contains transforms of all input fastas*/
	char **revTransform = bwt(Rm,seqCount,input.length);
	//printBwt(m,transform,seqCount);
	struct index *FMidx = calculateInterval(transform,input.length,seqCount,revTransform);
	//printInt(FMidx,seqCount,input.length);
	intervalToFile(FMidx,seqCount,m,transform,input);
    	deleteSuffixArray(m,seqCount,input.length);
    	deleteInputStruct(input,seqCount);
    	return 0;
}
