#pragma once
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include "matrix.h"

// dump some data to memory to be used by the generated code
template <class T>
int dumpToMem(vector<T>& data, const char* filePath, int numOfElements) {
  int fd;
  int result;
  //size_t size = data.size() * sizeof(T);
  size_t size = numOfElements * sizeof(T);

//  cout << "data size: " << data.size() << "\n";
//  cout << "size: " << size << "\n";

  fd = open(filePath, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
  if(fd == -1) {
    perror("Error opening file for writing");
    exit(EXIT_FAILURE);
  }

  result = lseek(fd, size-1, SEEK_SET);
  if(result == -1) {
    close(fd);
    perror("Error calling lseek() to 'stretch' the file");
    printf("%s\n",strerror(errno));
    exit(EXIT_FAILURE);
  }

  result = write(fd, "", 1);
  if(result != 1) {
    close(fd);
    perror("Error writing last byte of the file");
    exit(EXIT_FAILURE);
  }

  T* map = static_cast<T*>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(map == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }

  // TODO: vectorize
  for(int i = 0; i < numOfElements; i++) 
    map[i] = data[i];

  if(munmap(map, size) == -1)
    perror("Error un-mmapping the file");

  
  fd = open(filePath, O_RDWR);
  if(fd == -1) {
    perror("Error opening file for reading");
    exit(EXIT_FAILURE);
  }

  T* map_r = static_cast<T*>(mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if(map_r == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }

/*  cout << filePath << ":\n";
  for(int i = 0; i < numOfElements; i++) 
    cout << i << ": " << map[i] << "\n";
  cout << "\n";*/

  if(munmap(map, size) == -1)
    perror("Error un-mmapping the file");

  close(fd);

  return 0;
}

// This is a temporary function to read the output of calculateWorkload.py
void readWorkloadPartition(vector<vector<int>>& partitionList) {
  string fName("partition_932.csv");
  ifstream file(fName, ifstream::in);
  string line;

  // Iterate through each line and split the content using delimeter
  while(getline(file, line)){
    vector<int> vec;
    stringstream values(line);
    string val;
    while(getline(values,val,',')) {
      size_t found=val.find("\r");
      if (found!=string::npos)
        val.erase(found,1);

      if(val.empty())
        break;

      int value = 0;
      stringstream ss(val);
      ss >> value;
      vec.push_back(value);
    }
    partitionList.push_back(vec);
   /* cout << "level:\n";
    for(auto& v : vec)
      cout << v << ", ";
    cout << "\n";*/
  }

  file.close();
}

void closeParanthesis(string fileName, int prevPartCounter, int sigType) {
  std::fstream str(fileName + "/calculate" + to_string(prevPartCounter) + ".c");
  str.seekp (0, ios::end);
  str.write("\n}",2);  // put the paranthesis at the end

  if(sigType == 2) {            // update the signature
    str.seekp (0, ios::beg);
    str.write("//",2);  // comment out the signature type 1

    string firstLine;
    std::getline(str,firstLine);
    str.seekp (3+firstLine.size(), ios::beg);
    str.write("  ",2); // enable the signature type 2
  }

  str.close();
 }

void removeParanthesis(string fileName, int prevPartCounter) {
  std::fstream str(fileName + "/calculate" + to_string(prevPartCounter) + ".c");
  str.seekp (-3, ios::end);
  str.write(" ",1); // remove the closing paranthesis at the end
  str.close();
}

void printVect(vector<int>& vect) {
  for(auto& item: vect)
  std::cout << item << ", ";
  std::cout << "\n";
}

// VERIFICATION OF XREF
void verifyX(int rows, vector<pair<vector<int>, vector<int>>>& allParents, vector<vector<double>>& allRowValues, vector<double>& b, vector<double>& x) {
  vector<double> xVerify(rows, 0.000000);
  for(int i = 0; i < rows; i++) {
    vector<int>& parents = allParents[i].first;
    vector<double>& rowValues = allRowValues[i];

    double sum = 0.0;
    for(int j = 0 ; j < parents.size() ; j++)
      sum += rowValues[j] * x[parents[j]];
    xVerify[i]= (b[i] - sum)/rowValues.back();
//    printf("b[%d]: %.6f,  xVerify[%d]: %.6f\n",i, b[i], i, xVerify[i]);
  }

   int errCnt = 0;
 // printf("printing errors\n");
   for (int i = 0; i < rows ; i++)
     if(fabs(1.0000000 - xVerify[i]) >= 1e-2) {
       errCnt++;
  //     printf("x[%d]: %.5f\n",i,xVerify[i]);
     }

   printf("errCnt:%d\n", errCnt); 
}

