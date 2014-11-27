using namespace std;

// Include Section
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <fstream>


int* matrix;

inline int getIndex(unsigned int x, unsigned int y, unsigned int width){
	return y*width+x;
}

 void readInputFile(std::string inputFileName, int pos1[2], int pos2[2], 
		     std::string* outputString1, std::string* outputString2){
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



// Cost function gave in the Project Assignment
short cost(int x) {
	int i, n_iter = 20;
	double dcost = 0;
	for(i = 0; i < n_iter; i++)
		dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
	return (short) (dcost / n_iter + 0.1);
}

int lcs(string stringA, string stringB, int positionA, int positionB){
	int maxi = 0;
	int unsigned lengthA = stringA.length();
	int unsigned lengthB = stringB.length();
	for(unsigned int i=0; i<lengthB; i++){
	  #pragma omp parallel for shared(matrix)
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
		#pragma omp parallel for shared(matrix)
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

// Function that calculates the path for the lcs in the matrix
string runMatrix(string stringA, string stringB, int positionA, int positionB){
	// Local variables
	// return string
	string ret;
	// copy both positions to be able to change them saving the previous result
	int posA = positionA-1;
	int posB = positionB-1;

	// body
	// while the positions are in range
	while(posA >= 0 && posB >= 0){
		// if the chars match decrement both positions and add the char to the ret string
		if(stringA[posA] == stringB[posB]){
			ret = stringA[posA] + ret;
			posA--;
			posB--;
			// if the chars don't match
		} else{
			// if this is the position 0,0 of the matrix stop the loop
			if(posA==0 && posB==0)
				break;
			// if the A position is in the limit of the matrix decrement positionB
			else if(posA==0)
				posB--;
			// if the B position is in the limit of the matrix decrement positionA
			else if(posB==0)
				posA--;
			// if the A-1,B position value is larger then the A,B-1 position use that path
			else if(matrix[getIndex(posA-1,posB,positionA)] > matrix[getIndex(posA,posB-1,positionA)])
				posA--;
			// if the A,B-1 position value is larger then the A-1,B position use that path
			else if(matrix[getIndex(posA-1,posB,positionA)] <= matrix[getIndex(posA,posB-1,positionA)])
				posB--;
			// if this return is reached something went terribly wrong
			else return "sou um totó";
		}
	} // ends while loop
	// return ret string
	return ret;
}
/* MAIN: PROCESS PARAMETERS */
int main(int argc, char *argv[]) {
	omp_set_num_threads(8);
	// Local Variables
	int lengthA = 0;
	int lengthB = 0;
	string final = "";
	string filename = argv[1];
	string stringA;
	string stringB;

	// Create stream to read the file content
	std::ifstream infile(filename.c_str());

	// Scan input for lengths
	infile >> lengthA >> lengthB;


	// Scan input for the strings
	// 1st get line to read the \n
	getline(infile, stringA);
	getline(infile,stringA);
	getline(infile,stringB);

	matrix = (int*)malloc((lengthA)*(lengthB)*sizeof(int));
	std::cout << lcs(stringA, stringB, lengthA,lengthB) << endl;
	std::cout << runMatrix(stringA, stringB, lengthA,lengthB) << endl;
}
