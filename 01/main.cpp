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
    int nodesNum;
    int weightsSum;
    int edgesNum;
    int **matrix;
    ~GraphStruct() {
        for (int i = 0; i < nodesNum; i++) {
            delete [] matrix[i];
        }
        delete [] matrix;
    }
};

/**
 * @brief Structure represents the linked list node for edge 
 * 
 */
struct Edge {
    // edge data
    int weight;
    int n1;
    int n2;
    int isUsed; // 1 or 0 if it's already used
};


/**
 * @brief Structure represents the search result
 * 
 */
struct ResultNode {
    ResultNode* next;

    // data
    GraphStruct* graph;
    int* cNodes;
    int weight;
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
    for (int i = 0; i < graph->nodesNum; i++) {
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
    graph->nodesNum = nNodes;

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
 * @brief Function copies given graph
 * 
 * @param graph given graph
 * @return GraphStruct* new identical graph
 */
GraphStruct* copyGraph(GraphStruct* graph) {
    GraphStruct* newGraph = new GraphStruct;
    newGraph->matrix = new int*[graph->nodesNum];
    
    for (int i = 0; i < graph->nodesNum; i++) {
        newGraph->matrix[i] = new int[graph->nodesNum];
        for (int j = 0; j < graph->nodesNum; j++) {
            newGraph->matrix[i][j] = graph->matrix[i][j];
        }
    }

    newGraph->edgesNum = graph->edgesNum;
    newGraph->nodesNum = graph->nodesNum;
    newGraph->weightsSum = graph->weightsSum;
    
    return newGraph;
}

/**
 * @brief Function sums the edges weights
 * 
 * @param graph 
 * @return int 
 */
int sumEdgesWeights(GraphStruct* graph){
    int edgeSum = 0;
    for (int i = 0; i < graph->nodesNum; i++) {
        for (int j = i; j < graph->nodesNum; j++) {
            edgeSum += graph->matrix[i][j];
        }
    }
    return edgeSum;
}

/**
 * @brief Create array of sorted edges
 * 
 * @param graph graph with edges to be sorted
 * @return Edge** sorted array of edges
 */
Edge** createSorEdgesLL(GraphStruct* graph) {
    Edge** edges = new Edge*[graph->edgesNum];
    int edgeID = 0;

    // add the edges
    for (int i = 0; i < graph->nodesNum; i++) {
        for (int j = i; j < graph->nodesNum; j++) {
            edges[edgeID] = new Edge;
            edges[edgeID]->isUsed = 0;
            edges[edgeID]->weight = graph->matrix[i][j];
            edges[edgeID]->n1 = i;
            edges[edgeID]->n2 = j;
        }
    }

    // sort the edges - buble sort descending
    Edge * tmpEdge = edges[0];
    for (int i = 0; i < graph->nodesNum; i++) {
        tmpEdge = edges[i];
        for (int j = i+1; j < graph->nodesNum; j++) {
            if (edges[j]->weight > edges[i]->weight) {
                edges[i] = edges[j];
                edges[j] = tmpEdge;
            }
        }
    }

    return edges;
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
    for (int i = 0; i < graph->nodesNum; i++) {
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
    int* visNodes = new int[graph->nodesNum];
    
    // sets zeros to all nodes
    for (int i = 0; i < graph->nodesNum; i++) {
        visNodes[i] = 0;
    } 

    // search graph - start with 0
    srchdfs(graph, visNodes, 0);

    // test the visited nodes
    for (int i = 0; i < graph->nodesNum; i++) {
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
    for ( int i = 0; i  < graph->nodesNum; i++ ) {
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
    int* visNodes = new int[graph->nodesNum];
    for (int i = 0; i < graph->nodesNum; i++) {
        // set the color to undefined
        visNodes[i] = -1;
    }

    if (colorGraph(graph, visNodes, 0, 0)) {
        for (int i = 0; i < graph->nodesNum; i++) {
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
 * @brief Function adds new result to top od the linked list
 * 
 * @param results results so far, linked list
 * @param graph graph with result
 * @param cNodes colored nodes of given graph1
 * @return ResultNode* updated results
 */
ResultNode* addResult(ResultNode* results, GraphStruct* graph, int* cNodes) {
    ResultNode* newResult = new ResultNode;
    newResult->cNodes = cNodes;
    newResult->graph = graph;
    newResult->next = results;

    return newResult;
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
    graph->nodesNum = num;
    graph->matrix = new int*[num];
    for (int i = 0; i < num; i++) {
        graph->matrix[i] = new int[num];
    }

    int weightsSum = 0;
    int edgesNum = 0;
    
    // read the rest of the file
    for (int i = 0; i < num; i++) {
        getline(graphFile, line);
        istringstream iss(line);

        for (int j = 0; j < num; j++) {
            iss >> graph->matrix[i][j];
            if (graph->matrix[i][j] != 0){
                edgesNum++;
                weightsSum+=graph->matrix[i][j];
            }
        }
    }

    // fix the total weight and total edges num
    graph->edgesNum = edgesNum/2;
    graph->weightsSum = weightsSum/2;

    return graph;
}



ResultNode* getMaxBiparSubgraph(GraphStruct* graph) {
    // create all neccesary components
    ResultNode* results = NULL;
    int* cNodes = new int[graph->nodesNum];
    GraphStruct* subgraph = createGraph(graph->nodesNum);
    Edge** edges = createSorEdgesLL(graph);
    
    // begin the search

    return results;
}

/**
 * @brief Main funciton from witch the program runs
 * 
 * @param argc Number of input data
 * @param argv Input data
 * @return int 
 */
int main(int argc, char *argv[]) {
    // test if input is correct
    if (argc != 2) {
        cout << "Incorrect input, just name of graph is needed" << "\n";
        return 1;
    }

    string graphName = argv[1];
    GraphStruct* graph = loadGraph(graphName);

    if (graph == NULL) {
        return 0;
    }

    cout << graph->nodesNum << endl;

    cout << "weights " << graph->weightsSum << " edgs:" << graph->edgesNum << endl;

    // test if graph can be used
    if (!isConnected(graph)) {
        cout << "Graph is not connected" << endl;
        return 0;
    }
  
    // end if the graph itself is biparted
    if (isBiparted(graph)) {
        cout << "Graph is biparted" << endl;
    }

    ResultNode* results = getMaxBiparSubgraph(graph);

    delete graph;

    return 0;
}