//Get sequences as input into char array
//maybe use an integer array?
//find the commonalities between the arrays
//this part can be parallelized
//there must be communication between processes
//find all the potential primers
//combine all the potential primers

//how do you want to represent the sequences?
//2D integer array
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"


char intToCharNclt(int ncltInt);
int charToIntNclt(char ncltChar);

int main(int argc, char *argv[])
{
	
	int comm_sz;
	int my_rank;
	char hostname[128];
	size_t len = 126;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	gethostname(hostname,len);
	
	
	int MIN_PRIMER_SIZE = 2;
	int sequenceLen = 10;
	int numOfSeq = 2;
	int localseqsize = (sequenceLen) / comm_sz;
	int sequences [numOfSeq*sequenceLen];
	int localsequence [localseqsize * numOfSeq];
	int potentialPrIndex [sequenceLen];
	int localpotentialPrIndex [localseqsize];
	int primerList[sequenceLen][2];
	int potPrCt = 0;
	int glbPotPrCt = 0;
	int testSeq = -1;
	int i,j, k;
	double startT, endT;
	
	
	if (my_rank == 0){
		srand(time(NULL));
		for (i = 0; i < numOfSeq*sequenceLen; i++) {
			sequences[i] = rand() % 2;
		}
		
	}
	
	int totalSize = localseqsize*numOfSeq;		
	MPI_Scatter(sequences, totalSize, MPI_INT, localsequence, totalSize, MPI_INT,0, MPI_COMM_WORLD);

	//finding commonalities in all sequences
	for(i = 0; i < localseqsize; i++) {
		for(j = 0; j < numOfSeq; j++) {
			if (j%MIN_PRIMER_SIZE == 0)
				j=0;
			if (j == 0) {
				testSeq = localsequence[(my_rank*localseqsize) + i*j];
			}
			if (testSeq != localsequence[(my_rank*localseqsize) + i*j]) {
				break;
			}
			else if (j == numOfSeq-1) {
				localpotentialPrIndex[potPrCt] = my_rank * localseqsize + i;
				potPrCt++;
			}
		}
	}
	
	if (my_rank == 0){		
		int count = 0;
		startT = MPI_Wtime();
		for (i = 0; i < potPrCt; i++) {
			potentialPrIndex[count] = localpotentialPrIndex[i];
			count++;
		}
		for (i = 1; i < comm_sz; i++) {
			int localsz;
			MPI_Recv(&localsz, 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int localPrInd[localsz];
			MPI_Recv(localPrInd, localsz, MPI_INT, i,i+comm_sz, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (j = 0; j < localsz; j++) {
				potentialPrIndex[count] = localPrInd[j];
				count++;
			}
		}
	}
	else {
		MPI_Send(&potPrCt, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD);
		MPI_Send(localpotentialPrIndex, potPrCt, MPI_INT, 0, my_rank+comm_sz, MPI_COMM_WORLD);
	}
	MPI_Reduce(&potPrCt, &glbPotPrCt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank == 0) {
		for (i = 0; i < glbPotPrCt; i++) {
			printf("%d\n", potentialPrIndex[i]);
		}
		
		//Finds if all the potential nctds are next to each other
		int adjacent = 0;
		int start = 0;
		int primerList [sequenceLen][2];
		int primerCt = 0;
		int potPrNtd;
		
		for(i = 0; i < glbPotPrCt-1; i++) {
			potPrNtd = potentialPrIndex[i];
			if (adjacent == 0) {
				start = potPrNtd;
			}
			if (potPrNtd + 1 == potentialPrIndex[i+1]){
				adjacent++;
			}
			else {
				adjacent = 0;
			}

			if (adjacent == MIN_PRIMER_SIZE-1) {
				primerList[primerCt][0] = start;
				primerList[primerCt][1] = potentialPrIndex[i+1];
				adjacent = 0;
				start = i;
				primerCt++;
			}
		}
		
		printf("Primer Candidates: \n");
		for (i = 0; i < primerCt; i++) {
			printf("%d-%d\n",primerList[i][0],primerList[i][1]);
		}
		
		for (i = 0; i < primerCt; i++) {
			int start1 = primerList[i][0];
			int end1 = primerList[i][1];
			
			for (j=start1; j<=end1; j++) {
				char charNclt = intToCharNclt(sequences[j]);
				printf("%c",charNclt);
			}
			printf(" index:%d-%d\n",start1,end1);
		}
		endT = MPI_Wtime();
		printf("Time elapsed: %f", endT-startT);
	}
	
	MPI_Finalize();
	return 0;
}

char intToCharNclt(int ncltInt) {
	char charNcltList[20] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q',
		'R','S','T','V','W','Y'};
	if (ncltInt < 0 || ncltInt > 19){
		return '*';
	}
	char charNclt = charNcltList[ncltInt];
	return charNclt;
}

int charToIntNclt(char ncltChar) {
	
	char charNcltList[20] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q',
		'R','S','T','V','W','Y'};
	int i;
	for (i = 0; i < 20; i++) {
		if (ncltChar == charNcltList[i]){
			return i;
		}
	}
	return '*';
}
	
//upstream primers
//downstream primers
