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
struct GraphStruct {
    int gNum;
    int **matrix;
    ~GraphStruct() {
        for (int i = 0; i < gNum; i++) {
            delete [] matrix[i];
        }
        delete [] matrix;
    }
};

/**
 * @brief Structure represents the linked list node for edge 
 * 
 */
struct EdgeNode {
    EdgeNode* prevEdge;
    EdgeNode* nextEdge;

    // edge data
    int weight;
    int n1;
    int n2;

    ~EdgeNode() {
        if (prevEdge != NULL) delete prevEdge;
        if (nextEdge != NULL) delete nextEdge;
    }
};


/**
 * @brief Function checks if the node can be colored by given color
 * 
 * @param graph given graph with edges
 * @param cNodes colored nodes
 * @param cNode current node
 * @param color color to be used
 * @return true if node can be colored
 * @return false otherwise
 */
bool canBeColored(GraphStruct* graph, int* cNodes, int cNode, int color) {
    for (int i = 0; i < graph->gNum; i++) {
        if (graph->matrix[cNode][i] != 0 && cNodes[i] == color) return false;
    }
    return true;
}

/**
 * @brief Create a Graph object
 * 
 * @param nNodes number of nodes
 * @return GraphStruct* object with graph represantion
 */
GraphStruct* createGraph(int nNodes) {
    // if bad dimention return null pointer
    if (nNodes <= 0) return NULL;

    // create new graph
    GraphStruct* graph = new GraphStruct;

    // store the number of nodes
    graph->gNum = nNodes;

    // create new graph
    int** arr = new int*[nNodes];
    for (int i = 0; i < nNodes; i++) {
        arr[i] = new int[nNodes];
    }
    graph->matrix = arr;

    // return the new graph
    return graph;
}

/**
 * @brief Function sums the edges weights
 * 
 * @param graph 
 * @return int 
 */
int sumEdgesWeights(GraphStruct* graph){
    int edgeSum = 0;
    for (int i = 0; i < graph->gNum; i++) {
        for (int j = i; j < graph->gNum; j++) {
            edgeSum += graph->matrix[i][j];
        }
    }
    return edgeSum;
}

EdgeNode* createSorEdgesLL(GraphStruct* graph) {
    // todo
}

/**
 * @brief Function recursively search all available nodes and sotres
 *        them to given visNodes
 * 
 * @param graph representation of searched graph
 * @param visNodes visited nodes
 * @param curNode current node
 */
void srchdfs(GraphStruct* graph, int* visNodes, int curNode) {
    visNodes[curNode] = 1;
    for (int i = 0; i < graph->gNum; i++) {
        if ( visNodes[i] == 0 && graph->matrix[curNode][i] != 0) {
            srchdfs(graph, visNodes, i);
        }
    }
}

/**
 * @brief Function tests the Connectivity of given graph
 * 
 * @return true if Connectivity is satisfied
 * @return false otherwise
 */
bool isConnected(GraphStruct* graph) {
    // create graph unvisited nodes
    int* visNodes = new int[graph->gNum];
    
    // sets zeros to all nodes
    for (int i = 0; i < graph->gNum; i++) {
        visNodes[i] = 0;
    } 

    // search graph - start with 0
    srchdfs(graph, visNodes, 0);

    // test the visited nodes
    for (int i = 0; i < graph->gNum; i++) {
        // if even one node hasn't been visited  return false
        if (visNodes[i] == 0) {
            delete [] visNodes;
            return false;
        }
    }

    delete[] visNodes;

    // all nodes visited, retun true
    return true;
}

/**
 * @brief This function colors given graph by two colors
 * 
 * @param graph given graph
 * @param cNodes nodes to be colored
 * @param cNode current node id
 * @param color color to be used on node
 * @return true if the coloring was successful
 * @return false otherwise
 */
bool colorGraph(GraphStruct* graph, int* cNodes, int cNode, int color) {
    // check if the color can be used
    if (canBeColored(graph, cNodes, cNode, color) == false) return false;

    // color the node
    cNodes[cNode] = color;

    int nextColor = color == 0 ? 1 : 0;

    // continue the coloring
    for ( int i = 0; i  < graph->gNum; i++ ) {
        if (graph->matrix[cNode][i] != 0 && cNodes[i] == -1) {
            if (colorGraph(graph, cNodes, i, nextColor) == false) return false;
        }
    }

    return true;
}

/**
 * @brief Function tests the Bipartity of given graph
 * 
 * @return true if the Bipartity is satisfied
 * @return false otherwise
 */
bool isBiparted(GraphStruct* graph) {
    // create graph unvisited nodes
    int* visNodes = new int[graph->gNum];
    for (int i = 0; i < graph->gNum; i++) {
        // set the color to undefined
        visNodes[i] = -1;
    }

    if (colorGraph(graph, visNodes, 0, 0)) {
        for (int i = 0; i < graph->gNum; i++) {
            if (visNodes[i] == -1) {
                delete[] visNodes;
                return false;
            }
        }
    }
    else {
        delete[] visNodes;
        return false;
    }
    
    delete[] visNodes;

    return true;
}

/**
 * @brief Function loads the file with given name
 *        in the graph_mbp folder and stores its
 *        data to the array
 * 
 * @param fileName File name with graph prescription
 */
GraphStruct* loadGraph(string graphName) {
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
    GraphStruct * graph = new GraphStruct;
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
    GraphStruct* graph = loadGraph(graphName);

    cout << graph->gNum << endl;

    for (int i = 0; i < graph->gNum; i++){
        for (int j = 0; j < graph->gNum; j++) {
            cout << graph->matrix[i][j] << "|";
        }
        cout << endl;
    }

    if (isConnected(graph)) {
        cout << "Graph is connected" << endl;
    }
    else {
        cout << "Graph is not connected" << endl;
    }

    if (isBiparted(graph)) {
        cout << "Graph is biparted" << endl;
    }
    else {
        cout << "Graph is not biparted" << endl;
    }

    delete graph;

    return 0;
}