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
#include <sstream>

using namespace std;

/**
 * @brief Structure represents the graph
 * 
 */
struct graphStruct
{
    int gNum;
    int **matrix;
    ~graphStruct() {
        for (int i = 0; i < gNum; i++) {
            delete [] matrix[i];
        }
        delete [] matrix;
    }
};


/**
 * @brief Function tests the Connectivity of given graph
 * 
 * @return true if Connectivity is satisfied
 * @return false otherwise
 */
bool testConnectivity() {
    return false;
}

/**
 * @brief Function tests the Bipartity of given graph
 * 
 * @return true if the Bipartity is satisfied
 * @return false otherwise
 */
bool testBipartity() {
    //todo color the graph
    return false;
}

/**
 * @brief Function loads the file with given name
 *        in the graph_mbp folder and stores its
 *        data to the array
 * 
 * @param fileName File name with graph prescription
 */
graphStruct* loadGraph(string graphName) {
    // load the file
    string line;
    ifstream graphFile ("graf_mbp/" + graphName);
    
    if (!graphFile.is_open()) {
        cout << "Unable to open graf_mbp/" + graphName + "\n";
        return NULL;
    }

    // read the first line
    int num;
    if (getline (graphFile, line) ) {
        istringstream iss(line);

        iss >> num;
    }
    else cout << "Unable to read the first line.\n";

    cout << "your input is " << num << "\n";
    
    // create edge metrix
    graphStruct * graph = new graphStruct;
    graph->gNum = num;
    graph->matrix = new int*[num];
    for (int i = 0; i < num; i++) {
        graph->matrix[i] = new int[num];
    }

    // read the rest of the file
    for (int i = 0; i < num; i++) {
        getline(graphFile, line);
        istringstream iss(line);

        for (int j = 0; j < num; j++) {
            iss >> graph->matrix[i][j];    
        }
    }

    return graph;
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
    // test if input is correct
    if (argc != 2) {
        cout << "Incorrect input, just name of graph is needed" << "\n";
        return 1;
    }

    string graphName = argv[1];
    graphStruct* graph = loadGraph(graphName);

    cout << graph->gNum << endl;

    for (int i = 0; i < graph->gNum; i++){
        for (int j = 0; j < graph->gNum; j++) {
            cout << graph->matrix[i][j] << "|";
        }
        cout << endl;
    }

    delete graph;

    return 0;
}