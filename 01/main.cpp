/**
 * @file main.cpp
 * @author David Omrai (omraidav@fit.cvut.com)
 * @brief First ni-pdp project sequence algorithm
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

/**
 * @brief Function loads the file with given name
 *        in the graph_mbp folder and stores its
 *        data to the array
 * 
 * @param fileName File name with graph prescription
 */
void loadGraph(string fileName) {
    // load the file
    string line;
    ifstream graphFile ("graf_mbp/" + fileName);

    if (!graphFile.is_open()) {
        cout << "Unable to open graf_mbp/" + fileName + "\n";
        return;
    }

    // load the first line
    if (getline (graphFile, line) ) {
        cout << line << "\n";
    }
    else cout << "Unable to read the first line.\n";

    // read the rest of the file
    
}

/**
 * @brief Main funciton from witch the program runs
 * 
 * @param argc Number of input data
 * @param argv Input data
 * @return int 
 */
int main(int argc, char *argv[]) {
    // todo
    cout << "Hello world!";
    return 0;
}