using namespace std;

// Include Section
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <fstream>
#include <mpi.h>

#define DEBUG 0

inline int getIndex(unsigned int x, unsigned int y, unsigned int width){
	return y*width+x;
}

// Cost function gave in the Project Assignment
short cost(int x) {
	int i, n_iter = 20;
	double dcost = 0;
	for(i = 0; i < n_iter; i++)
		dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
	return (short) (dcost / n_iter + 0.1);
}

/*
int lcs(string stringA, string stringB, int positionA, int positionB){
	int maxi = 0;
	int unsigned lengthA = stringA.length();
	int unsigned lengthB = stringB.length();
	for(unsigned int i=0; i<lengthB; i++){
	  //#pragma omp parallel for shared(matrix)
	  for(int unsigned j=0; j<min((unsigned int)i+1,lengthA); j++){
			unsigned int x = j;
			unsigned int y = i - j;
			if(stringA[x] == stringB[y]){
				if(x == 0 || y == 0)
					matrix[getIndex(x,y,positionA)] = 1;
				else matrix[getIndex(x,y,positionA)] = matrix[getIndex(x-1,y-1,positionA)] + cost(1);
			} else if(x == 0 && (y == 0)) matrix[getIndex(x,y,positionA)] = 0;
			else if(x==0)
				matrix[getIndex(x,y,positionA)] = matrix[getIndex(x,y-1,positionA)];
			else if(y==0)
				matrix[getIndex(x,y,positionA)] = matrix[getIndex(x-1,y,positionA)];
			else matrix[getIndex(x,y,positionA)] = std::max(matrix[getIndex(x-1,y,positionA)], matrix[getIndex(x,y-1,positionA)]);

			if(matrix[getIndex(x,y,positionA)]>maxi) maxi = matrix[getIndex(x,y,positionA)];
		}
	}

	for(int i=1; i<lengthA; i++){
		unsigned int lim = lengthA - i;
		lim = min(lim, lengthB);
		//#pragma omp parallel for shared(matrix)
		for(int j=0; j < lim; j++){
			int x = i + j;
			int y = lengthB - j - 1;
			if(stringA[x] == stringB[y]){
				if(x == 0 || (y == 0))
					matrix[getIndex(x,y,positionA)] = 1;
				else matrix[getIndex(x,y,positionA)] = matrix[getIndex(x-1,y-1,positionA)] + cost(1);
			} else if(x == 0 && y == 0) matrix[getIndex(x,y,positionA)] = 0;
			else if(x==0)
				matrix[getIndex(x,y,positionA)] = matrix[getIndex(x,y-1,positionA)];
			else if(y==0)
				matrix[getIndex(x,y,positionA)] = matrix[getIndex(x-1,y,positionA)];
			else matrix[getIndex(x,y,positionA)] = std::max(matrix[getIndex(x-1,y,positionA)], matrix[getIndex(x,y-1,positionA)]);

			if(matrix[getIndex(x,y,positionA)]>maxi) maxi = matrix[getIndex(x,y,positionA)];
		}
	}
	return maxi;
}
*/

int* computeMatrixBlock (string stringA, string stringB, int* dependencyList){
	int unsigned lengthA = stringA.length();
	int unsigned lengthB = stringB.length();
	int* matrix = (int*)malloc((lengthA)*(lengthB)*sizeof(int));

	for(unsigned int i=0; i<lengthB; i++){
	  //#pragma omp parallel for shared(matrix)
		for(int unsigned j=0; j<min((unsigned int)i+1,lengthA); j++){
			unsigned int x = j;
			unsigned int y = i - j;

			if(stringA[x] == stringB[y]){
				if(x == 0 && (y == 0)){
					matrix[getIndex(x,y,lengthA)] = dependencyList[0]+1;
				} else if(x==0){
					matrix[getIndex(x,y,lengthA)] = dependencyList[lengthA+y]+1;
				}else if(y==0){
					matrix[getIndex(x,y,lengthA)] = dependencyList[x]+1;
				}else {
					matrix[getIndex(x,y,lengthA)] = matrix[getIndex(x-1,y-1,lengthA)]+1;
				}
			}else{
				if(x == 0 && (y == 0)){
					matrix[getIndex(x,y,lengthA)] = max(dependencyList[1],dependencyList[lengthA+1]);
				} else if(x==0){
					matrix[getIndex(x,y,lengthA)] = max(matrix[getIndex(x,y-1,lengthA)], dependencyList[lengthA+y]);
				}else if(y==0){
					matrix[getIndex(x,y,lengthA)] = max(matrix[getIndex(x-1,y,lengthA)], dependencyList[x]);
				}else {
					matrix[getIndex(x,y,lengthA)] = max(matrix[getIndex(x-1,y,lengthA)], matrix[getIndex(x,y-1,lengthA)]);
				}
			}
		}
	}

	for(int i=1; i<lengthA; i++){
		unsigned int lim = lengthA - i;
		lim = min(lim, lengthB);
		//#pragma omp parallel for shared(matrix)
		for(int j=0; j < lim; j++){
			int x = i + j;
			int y = lengthB - j - 1;

			if(stringA[x] == stringB[y]){
				if(x == 0 && (y == 0)){
					matrix[getIndex(x,y,lengthA)] = dependencyList[0]+1;
				} else if(x==0){
					matrix[getIndex(x,y,lengthA)] = dependencyList[lengthA+y]+1;
				}else if(y==0){
					matrix[getIndex(x,y,lengthA)] = dependencyList[x]+1;
				}else {
					matrix[getIndex(x,y,lengthA)] = matrix[getIndex(x-1,y-1,lengthA)]+1;
				}
			}else{
				if(x == 0 && (y == 0)){
					matrix[getIndex(x,y,lengthA)] = max(dependencyList[1],dependencyList[lengthA+1]);
				} else if(x==0){
					matrix[getIndex(x,y,lengthA)] = max(matrix[getIndex(x,y-1,lengthA)], dependencyList[lengthA+y]);
				}else if(y==0){
					matrix[getIndex(x,y,lengthA)] = max(matrix[getIndex(x-1,y,lengthA)], dependencyList[x]);
				}else {
					matrix[getIndex(x,y,lengthA)] = max(matrix[getIndex(x-1,y,lengthA)], matrix[getIndex(x,y-1,lengthA)]);
				}
			}
		}
	}

	return matrix;
}

void printCurrentProcess(int* coords, int* initialPositions, int* finalPositions, int id){
	std::cout << "Process " << id << " starts running matrix block " << coords[0] << "," << coords[1]
	<< " with real values " << initialPositions[0] << "," << initialPositions[1] << " to "
	<< finalPositions[0] << "," << finalPositions[1] << "." << std::endl;
}

void readInputFile(std::string inputFileName, int pos1[2], int pos2[2], std::string* outputString1, std::string* outputString2){
	std::ifstream infile(inputFileName.c_str());
	char c;
	int counter = 0;
    //ignoring first line
	while(c!='\n'){
		infile.get(c);
	}

    //ignoring chars until initial position of first string
	while(counter < pos1[0]){
		infile.get(c);
		counter++;
	}
    //pushing back chars to string1
	for(int i=0; i < pos1[1] - pos1[0]; i++){
		infile.get(c);
		outputString1->push_back(c);
		infile.get(c);
	}
	counter=0;
    //ignoring chars until first position of second string
	while(counter < pos2[0]){
		infile.get(c);
		counter++;
	}
    //pushing back chars to string2
	for(int i=0; i < pos2[1] - pos2[0]; i++){
		infile.get(c);
		outputString2->push_back(c);
	}
	infile.close();
}

void getTargetDimension(int* targetCoords, int* dim, int lengthA, int lengthB, int* targetColumns, int* targetRows){
	int initialPositions[2];
	int finalPositions[2];
	initialPositions[0] = targetCoords[0]*(lengthA/dim[0]);
	initialPositions[1] = targetCoords[1]*(lengthB/dim[1]);
	finalPositions[0] = ((targetCoords[0]+1)*(lengthA/dim[0]))-1;
	finalPositions[1] = ((targetCoords[1]+1)*(lengthB/dim[1]))-1;
	if(finalPositions[0] > lengthA-1) finalPositions[0] = lengthA;
	if(finalPositions[1] > lengthB-1) finalPositions[1] = lengthB;
	*targetColumns = finalPositions[0] - initialPositions[0];
	*targetRows = finalPositions[1] - initialPositions[1];
	return;
}

void routineMPI(string filename, int n, int m){

	string stringA = "";
	string stringB = "";

	// Initialize MPI
	int id, num_procs;
	// write number of processes in num_procs
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	// write the current process identifier in id
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	MPI_Comm cart_comm;
	MPI_Status status;

	int initialPositions[2];
	int finalPositions[2];
	int dim[2];
	int periodic[2];
	int coords[2];

	memset(dim, 0, 2*sizeof(int));
	memset(periodic, 0, 2*sizeof(int));
	memset(coords, 0, 2*sizeof(int));

	MPI_Dims_create(num_procs, 2, dim);
	//Create Comunication
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periodic, 1, &cart_comm);

	if(DEBUG && id==0) std::cout << "Block matrix dimension: " << dim[0] << "," << dim[1] << std::endl;

	// Check for coordinates
	MPI_Cart_coords(cart_comm, id, 2, coords);
	// Calculate real coordinates
	initialPositions[0] = coords[0]*(n/dim[0]);
	initialPositions[1] = coords[1]*(m/dim[1]);
	finalPositions[0] = ((coords[0]+1)*(n/dim[0]))-1;
	finalPositions[1] = ((coords[1]+1)*(m/dim[1]))-1;
	if(finalPositions[0] > n-1) finalPositions[0] = n;
	if(finalPositions[1] > m-1) finalPositions[1] = m;

	// Initialize local values to each process
	int maxColumn = finalPositions[0] - initialPositions[0];
	int maxRow = finalPositions[1] - initialPositions[1];
	int dependencies = maxColumn + maxRow + 1;

	// Matrix declaration
	int* matrix;

	readInputFile(filename, initialPositions, finalPositions, &stringA, &stringB);

	if(coords[0]==0 && coords[1]==0){
		if(DEBUG) printCurrentProcess(coords, initialPositions, finalPositions, id);

		if(dim[0] > 1){
			int targetSideCoords[2];
			targetSideCoords[0] = 1;
			targetSideCoords[1] = 0;

			int targetSideId;
			int targetSideColumns, targetSideRows;

			getTargetDimension(targetSideCoords, dim, n, m, &targetSideColumns, &targetSideRows);

			int dependencyListSideSize = targetSideColumns+targetSideRows+1;
			int dependencyListSide[dependencyListSideSize];

		// set dependencies of Side block
			for(int i = 0; i<targetSideColumns+1; i++)
				dependencyListSide[i]=0;
			for(int i = targetSideColumns+1; i< dependencyListSideSize; i++)
				dependencyListSide[i]=matrix[getIndex(maxColumn-1,i-targetSideColumns+1, maxColumn)];

		// send to side
			MPI_Cart_rank(cart_comm, targetSideCoords, &targetSideId);
			if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetSideId << " with length " << dependencyListSideSize << " with tag " << targetSideId << std::endl << std::flush;
			MPI_Send(dependencyListSide, dependencyListSideSize, MPI_INT, targetSideId, targetSideId, MPI_COMM_WORLD);
		}

		if(dim[1] > 1){
			int targetBotCoords[2];
			targetBotCoords[0] = 0;
			targetBotCoords[1] = 1;

			int targetBotId;
			int targetBotColumns, targetBotRows;

			getTargetDimension(targetBotCoords, dim, n, m, &targetBotColumns, &targetBotRows);

			int dependencyListBotSize = targetBotColumns+targetBotRows+1;

			int dependencyListBot[dependencyListBotSize];

		// set dependencies of Bot block
			dependencyListBot[0]=0;
			for(int i = 1; i<targetBotColumns+1; i++)
				dependencyListBot[i]=matrix[getIndex(i-1,maxRow-1, maxColumn)];
			for(int i = targetBotColumns+1; i < dependencyListBotSize; i++)
				dependencyListBot[i] = 0;

		// send to bot
			MPI_Cart_rank(cart_comm, targetBotCoords, &targetBotId);
			MPI_Send(dependencyListBot, dependencyListBotSize, MPI_INT, targetBotId, targetBotId, MPI_COMM_WORLD);
			if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetBotId << " with length " << dependencyListBotSize << " with tag " << targetBotId << std::endl << std::flush;
		}
	} else{
		// is this matrix block in the first column or first row?
		if(coords[0]==0){ // matrix block in the first column
			int myDependencies = maxColumn + maxRow + 1;
			int dependencyList[myDependencies];

			int sourceCoords[2];
			sourceCoords[0] = coords[0];
			sourceCoords[1] = coords[1] - 1;

			if(DEBUG) printCurrentProcess(coords, initialPositions, finalPositions, id);

			int sourceId;
			MPI_Cart_rank(cart_comm, sourceCoords, &sourceId);

			if(DEBUG) std::cout << "Process " << id << " is waiting message from process " << sourceId << " with length " << myDependencies << " and tag " << id << std::endl << std::flush;
			MPI_Recv(dependencyList, myDependencies, MPI_INT, sourceId, id, MPI_COMM_WORLD, &status);
			if(DEBUG) std::cout << "Process " << id << " received message from process " << sourceId << std::endl << std::flush;

			//run function matrix(...) TODO

			// send to bot block
			if(coords[1] < dim[1]-1){
				int targetBotCoords[2];
				targetBotCoords[0] = coords[0];
				targetBotCoords[1] = coords[1] + 1;

				int targetBotId;
				int targetBotColumns, targetBotRows;

				getTargetDimension(targetBotCoords, dim, n, m, &targetBotColumns, &targetBotRows);

				int dependencyListBotSize = targetBotColumns+targetBotRows+1;

				int dependencyListBot[dependencyListBotSize];

				dependencyListBot[0] = 0;
				for(int i = 1; i<targetBotColumns+1; i++)
					dependencyListBot[i]=matrix[getIndex(i-1,maxRow-1, maxColumn)];
				for(int i = targetBotColumns+1; i < dependencyListBotSize; i++)
					dependencyListBot[i] = 0;

				if(DEBUG) std::cout << "Process " << id << " is trying to find " << targetBotCoords[0] << "," << targetBotCoords[1] << std::endl << std::flush;
				MPI_Cart_rank(cart_comm, targetBotCoords, &targetBotId);
				if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetBotId << " with length " << dependencyListBotSize << std::endl << std::flush;
				MPI_Send(dependencyListBot, dependencyListBotSize, MPI_INT, targetBotId, targetBotId, MPI_COMM_WORLD);
			}

			if(coords[0] < dim[0]-1){
			// create destinations links and send data
				int targetSideCoords[2];
				targetSideCoords[0] = coords[0] + 1;
				targetSideCoords[1] = coords[1];

				int targetSideId;
				int targetSideColumns, targetSideRows;

				getTargetDimension(targetSideCoords, dim, n, m, &targetSideColumns, &targetSideRows);

			//create dependency list
				int dependencyListSide[targetSideRows];

				for(int i = 0; i < targetSideRows; i++){
					dependencyListSide[i]=matrix[getIndex(maxColumn-1,i, maxColumn)];
				}

			// send to Side block
				if(DEBUG) std::cout << "Process " << id << " is trying to find " << targetSideCoords[0] << "," << targetSideCoords[1] << std::endl << std::flush;
				MPI_Cart_rank(cart_comm, targetSideCoords, &targetSideId);
				if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetSideId << " with length " << targetSideRows << std::endl << std::flush;
				MPI_Send(dependencyListSide, targetSideRows, MPI_INT, targetSideId, targetSideId, MPI_COMM_WORLD);
			}
		} else if(coords[1]==0){ //matrix block in the first row
			int myDependencies = maxColumn + maxRow + 1;
			int dependencyList[myDependencies];

			int sourceCoords[2];
			sourceCoords[0] = coords[0] - 1;
			sourceCoords[1] = coords[1];

			if(DEBUG) printCurrentProcess(coords, initialPositions, finalPositions, id);

			int sourceId;
			MPI_Cart_rank(cart_comm, sourceCoords, &sourceId);

			if(DEBUG) std::cout << "Process " << id << " is waiting message from process " << sourceId << " with length " << myDependencies << " and tag " << id << std::endl << std::flush;
			MPI_Recv(dependencyList, myDependencies, MPI_INT, sourceId, id, MPI_COMM_WORLD, &status);
			if(DEBUG) std::cout << "Process " << id << " received message from process " << sourceId << std::endl << std::flush;

			//run function matrix(...) TODO
			if(coords[1] < dim[1]-1){
			// create destinations links and send data
				int targetBotCoords[2];
				targetBotCoords[0] = coords[0];
				targetBotCoords[1] = coords[1] + 1;

				int targetBotId;
				int targetBotColumns, targetBotRows;

				getTargetDimension(targetBotCoords, dim, n, m, &targetBotColumns, &targetBotRows);

			//create both dependency lists
				int dependencyListBot[targetBotColumns + 1];

				dependencyListBot[0] = dependencyList[myDependencies-1];
				for(int i = 1; i<targetBotColumns+1; i++)
					dependencyListBot[i]=matrix[getIndex(i-1,maxRow-1, maxColumn)];

				if(DEBUG) std::cout << "Process " << id << " is trying to find " << targetBotCoords[0] << "," << targetBotCoords[1] << std::endl << std::flush;
				MPI_Cart_rank(cart_comm, targetBotCoords, &targetBotId);
				if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetBotId << " with length " << targetBotColumns + 1 << " with tag " << targetBotId << std::endl << std::flush;
				MPI_Send(dependencyListBot, targetBotColumns + 1, MPI_INT, targetBotId, targetBotId, MPI_COMM_WORLD);
			}

			if(coords[0] < dim[0]-1){
				int targetSideCoords[2];
				targetSideCoords[0] = coords[0] + 1;
				targetSideCoords[1] = coords[1];

				int targetSideId;
				int targetSideColumns, targetSideRows;

				getTargetDimension(targetSideCoords, dim, n, m, &targetSideColumns, &targetSideRows);

				int dependencyListSideSize = targetSideColumns+targetSideRows+1;
				int dependencyListSide[dependencyListSideSize];

				for(int i = 0; i < targetSideColumns+1; i++){
					dependencyListSide[i]=0;
				}
				for(int i = targetSideColumns+1; i < dependencyListSideSize; i++){
					dependencyListSide[i]=matrix[getIndex(maxColumn-1,i - (targetSideColumns+1), maxColumn)];
				}
				if(DEBUG) std::cout << "Process " << id << " is trying to find " << targetSideCoords[0] << "," << targetSideCoords[1] << std::endl << std::flush;
				MPI_Cart_rank(cart_comm, targetSideCoords, &targetSideId);
				if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetSideId << " with length " << dependencyListSideSize << std::endl << std::flush;
				MPI_Send(dependencyListSide, dependencyListSideSize, MPI_INT, targetSideId, targetSideId, MPI_COMM_WORLD);
			}

		} else{
			int myDependencies = maxColumn + maxRow + 1;
			int myDependencyList[myDependencies];

			int sourceCoordsTop[2];
			int sourceCoordsSide[2];

			sourceCoordsTop[0] = coords[0];
			sourceCoordsTop[1] = coords[1] - 1;

			sourceCoordsSide[0] = coords[0] - 1;
			sourceCoordsSide[1] = coords[1];

			if(DEBUG) printCurrentProcess(coords, initialPositions, finalPositions, id);

			int sourceIdTop, sourceIdSide;
			MPI_Cart_rank(cart_comm, sourceCoordsTop, &sourceIdTop);
			MPI_Cart_rank(cart_comm, sourceCoordsSide, &sourceIdSide);

			int myDependencyListSide[maxRow];
			int myDependencyListTop[maxColumn + 1];

			if(DEBUG) std::cout << "Process " << id << " is waiting message from process " << sourceIdSide << " with length " << maxRow << std::endl << std::flush;
			MPI_Recv(myDependencyListSide, maxRow, MPI_INT, sourceIdSide, id, MPI_COMM_WORLD, &status);
			if(DEBUG) std::cout << "Process " << id << " received message from process " << sourceIdSide << std::endl << std::flush;

			if(DEBUG) std::cout << "Process " << id << " is waiting message from process " << sourceIdTop << " with length " << maxColumn + 1 << std::endl << std::flush;
			MPI_Recv(myDependencyListTop, maxColumn + 1, MPI_INT, sourceIdTop, id, MPI_COMM_WORLD, &status);
			if(DEBUG) std::cout << "Process " << id << " received message from process " << sourceIdTop << std::endl << std::flush;

			for(int i = 0; i<maxColumn+1; i++){
				myDependencyList[i] = myDependencyListTop[i];
			}
			for(int i = maxColumn+1; i < myDependencies; i++){
				myDependencyList[i] = myDependencyListSide[i - (maxColumn+1)];
			}


			if(coords[0] < dim[0]-1){
				int targetSideCoords[2];
				targetSideCoords[0] = coords[0] + 1;
				targetSideCoords[1] = coords[1];

				int targetSideId;
				int targetSideColumns, targetSideRows;

				getTargetDimension(targetSideCoords, dim, n, m, &targetSideColumns, &targetSideRows);

				int dependencyListSideSize = targetSideColumns+targetSideRows+1;

				int dependencyListSide[targetSideRows];

				for(int i = 0; i < targetSideRows; i++){
					dependencyListSide[i]=matrix[getIndex(maxColumn-1,i, maxColumn)];
				}

				MPI_Cart_rank(cart_comm, targetSideCoords, &targetSideId);
				if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetSideId << " with length " << targetSideRows << std::endl << std::flush;
				MPI_Send(dependencyListSide, targetSideRows, MPI_INT, targetSideId, targetSideId, MPI_COMM_WORLD);
			}

			if(coords[1] < dim[1]-1){
				int targetBotCoords[2];
				targetBotCoords[0] = coords[0];
				targetBotCoords[1] = coords[1] + 1;

				int targetBotId;
				int targetBotColumns, targetBotRows;

				getTargetDimension(targetBotCoords, dim, n, m, &targetBotColumns, &targetBotRows);

				int dependencyListBotSize = targetBotColumns+targetBotRows+1;

				int dependencyListBot[targetBotColumns + 1];

				dependencyListBot[0] = myDependencyList[myDependencies];
				for(int i = 1; i < targetBotColumns+1; i++){
					dependencyListBot[i]=matrix[getIndex(i-1,maxRow-1, maxColumn)];
				}

				MPI_Cart_rank(cart_comm, targetBotCoords, &targetBotId);
				if(DEBUG) std::cout << "Process " << id << " sent message to process " << targetBotId << " with length " << targetBotColumns + 1 << std::endl  << std::flush;
				MPI_Send(dependencyListBot, targetBotColumns + 1, MPI_INT, targetBotId, targetBotId, MPI_COMM_WORLD);
			}

			if(coords[0]==dim[0]-1 && coords[1]==dim[1]-1)
				std::cout << "last block" << std::endl;
		}
	}
}

/*
void runLocalMatrix(string stringA, string stringB, int positionA, int positionB, MPI_Comm cart_comm, int coords[2], string* ret){
	int posA = positionA-1;
	int posB = positionB-1;
	int destinationCoords[2];
	int destinationId;
	int localMatrixCoords[2]; // traduzir local de um bloco para o outro
	int destMatrixCoords[2];

	while(posA >= 0 && posB >= 0){
		localMatrixCoords[0] = posA;
		localMatrixCoords[1] = posB;
			// if the chars match decrement both positions and add the char to the ret string
		if(stringA[posA] == stringB[posB]){
			*ret = stringA[posA] + *ret;
				if(posA == 0 && posB == 0){ //mandar para o bloco diagonal
					destinationCoords[0] = coords[0]-1;
					destinationCoords[1] = coords[1]-1;
					MPI_Cart_rank(cart_comm, destinationCoords, &destinationId);
						destMatrixCoords[0]=//maxTamanho de colunas da matriz destino
						destMatrixCoords[1]=//maxTamanho de colunas da matriz destino
						MPI_Send(destMatrixCoords, 2, MPI_INT, destinationId, destinationId, MPI_COMM_WORLD); //send matrix position
						//MPI_Send(ret, ret.length(), MPI_CHAR, destinationId, destinationId, MPI_COMM_WORLD); //send constructed string
					}else{
					if(posA == 0){ //mandar para o bloco da esquerda na diagonal
						destinationCoords[0] = coords[0];
						destinationCoords[1] = coords[1]-1;
						MPI_Cart_rank(cart_comm, destinationCoords, &destinationId);
						destMatrixCoords[0]=//maxTamanho de colunas da matriz destino
						destMatrixCoords[1]=localMatrixCoords[1]-1;
						MPI_Send(destMatrixCoords, 2, MPI_INT, destinationId, destinationId, MPI_COMM_WORLD); //send matrix position
						//MPI_Send(ret, ret.length(), MPI_CHAR, destinationId, destinationId, MPI_COMM_WORLD); //send constructed string
					}else if(posB == 0){ //mandar para o bloco de cima na diagonal
						destinationCoords[0] = coords[0]-1;
						destinationCoords[1] = coords[1];
						MPI_Cart_rank(cart_comm, destinationCoords, &destinationId);
						destMatrixCoords[0]=localMatrixCoords[0]-1;
						destMatrixCoords[1]=//maxTamanho de linhas da matriz destino
						MPI_Send(destMatrixCoords, 2 , MPI_INT, destinationId, destinationId, MPI_COMM_WORLD);
						//MPI_Send(ret,ret.length() , MPI_CHAR, destinationId, destinationId, MPI_COMM_WORLD); //send constructed string
					}else{
						posA--;
						posB--;
					}
				}
			} else{
				if(coords[0]==0 && coords[1]==0){	//se estiver no bloco canto superior esquerdo
					if(posA==0 && posB==0)
						break;
					else if(posA==0)
						posB--;
					else if(posB==0)
						posA--;
					else if(matrix[getIndex(posA-1,posB,positionA)] > matrix[getIndex(posA,posB-1,positionA)])
						posA--;
					else if(matrix[getIndex(posA-1,posB,positionA)] <= matrix[getIndex(posA,posB-1,positionA)])
						posB--;
				}else if(coords[0]==0 || coords[1]==0){ //se estiver ou na linha ou na coluna de blocos limite
					if(coords[0]==0){ //linha dos blocos mais a cima de todos
						if(posB==0)
							posA--;
						else if(matrix[getIndex(posA-1,posB,positionA)] > matrix[getIndex(posA,posB-1,positionA)])
							posA--;
						else if(matrix[getIndex(posA-1,posB,positionA)] <= matrix[getIndex(posA,posB-1,positionA)])
							posB--;
					}
					if(coords[1]==0){ //coluna dos blocos mais a esquerda de todos
						if(posA==0)
							posB--;
						else if(matrix[getIndex(posA-1,posB,positionA)] > matrix[getIndex(posA,posB-1,positionA)])
							posA--;
						else if(matrix[getIndex(posA-1,posB,positionA)] <= matrix[getIndex(posA,posB-1,positionA)])
							posB--;
					}
				}else{	//caso seja um bloco que esteja no interior da matriz de blocos
					if(posA==0){	//caso esteja na coluna limite esquerda da matriz local
						destinationCoords[0] = coords[0];
						destinationCoords[1] = coords[1]-1;
						MPI_Cart_rank(cart_comm, destinationCoords, &destinationId);
						destMatrixCoords[0]=//maxTamanho de colunas da matriz destino
						destMatrixCoords[1]=localMatrixCoords[1];
						MPI_Send(destMatrixCoords, 2, MPI_INT, destinationId, destinationId, MPI_COMM_WORLD); //send matrix position
						//MPI_Send(ret, ret.length(), MPI_CHAR, destinationId, destinationId, MPI_COMM_WORLD); //send constructed string
					}else if(matrix[getIndex(posA-1,posB,positionA)] > matrix[getIndex(posA,posB-1,positionA)])
					posA--;
					else if(matrix[getIndex(posA-1,posB,positionA)] <= matrix[getIndex(posA,posB-1,positionA)])
						posB--;
				}
			}
		}// ends while loop
	}


	string runMatrix(int dim[2],string stringA, string stringB,int positionA, int positionB){
	// Local variables
	// return string
		int id;
		string ret;
		MPI_Comm cart_comm;
		MPI_Status status;
		int coords[2];
		int buf[2];
		MPI_Cart_coords(cart_comm, id, 2, coords);
	// body
	// while the positions are in range
		if(coords[0] == dim[0]-1 && coords[1] == dim[1]-1){
			runLocalMatrix(stringA, stringB, positionA, positionB, cart_comm, coords, &ret);
		}else{
		MPI_Recv(buf,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status); //CHANGE TAG
		runLocalMatrix(stringA, stringB, positionA, positionB, cart_comm, coords, &ret);
	}
	// return ret string
	return ret;
}*/

/* MAIN: PROCESS PARAMETERS */
	int main(int argc, char *argv[]) {
	//omp_set_num_threads(8);

		string filename = argv[1];
		int n, m;

	// Create stream to read the file content
		std::ifstream infile(filename.c_str());
	// Scan input for lengths
		infile >> n >> m;
		infile.close();

		MPI_Init(&argc,&argv);
		routineMPI(filename, n, m);
		MPI_Finalize();
	//

	//std::cout << lcs(stringA, stringB, lengthA,lengthB) << endl;
	//std::cout << runMatrix(dim,stringA, stringB, lengthA,lengthB) << endl;
	}
