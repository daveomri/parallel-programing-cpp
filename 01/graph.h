/**
 * @file graph.h
 * @author David Omrai
 * @brief Functions for operations with graph
 * @version 0.1
 * @date 2022-03-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include <iostream>
#include <sstream>  
#include <string>
#include <fstream>

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
    ~ResultNode() {
        delete graph;
        delete [] cNodes;
    }
};

struct Results
{
    ResultNode* results;
};

// ----------------------------------------------------

/**
 * @brief Function loads the file with given name
 *        in the graph_mbp folder and stores its
 *        data to the array
 * 
 * @param fileName File name with graph prescription
 */
GraphStruct* loadGraph(std::string graphName) {
    // load the file
    std::string line;
    std::ifstream graphFile ("graf_mbp/" + graphName);
    
    if (!graphFile.is_open()) {
        std::cout << "Unable to open graf_mbp/" + graphName + "\n";
        return NULL;
    }

    // read the first line
    int num;
    if (getline (graphFile, line) ) {
        std::istringstream iss(line);

        iss >> num;
    }
    else std::cout << "Unable to read the first line.\n";
    
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
        std::getline(graphFile, line);
        std::istringstream iss(line);

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

// ----------------------------------------------------


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

    // store the data
    graph->nodesNum = nNodes;
    graph->weightsSum = 0;
    graph->edgesNum = 0;

    // create new graph
    int** arr = new int*[nNodes];
    for (int i = 0; i < nNodes; i++) {
        arr[i] = new int[nNodes];
        for (int j = 0; j < nNodes; j++) {
            arr[i][j] = 0; 
        }
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
 * @brief Function copies given vector
 * 
 * @param vec vector to be copies
 * @param len length of vector
 * @return int* copy of given vector
 */
int* copyVector(int* vec, int len) {
    int* copyVec = new int[len];

    for (int i = 0; i < len; i++) {
        copyVec[i] = vec[i];
    }

    return copyVec;
}

/**
 * @brief Function tests if the node can be colored with given color
 * 
 * @param graph given graph
 * @param cNodes array of colored nodes
 * @param cNode id of curren node
 * @param color color to be used on cNode
 * @return true if can be colored
 * @return false otherwise
 */
bool canColorNode(GraphStruct* graph, int* cNodes, int cNode, int color) {
    for (int i = 0; i < graph->nodesNum; i++) {
        if (graph->matrix[cNode][i] != 0 && cNodes[i] == color) return false;
    }
    return true;
}

/**
 * @brief Function adds edge to graph
 * 
 * @param graph given graph
 * @param edge edge to be added
 */
void addEdge(GraphStruct* graph, Edge* edge) {
    graph->edgesNum += 1;
    graph->weightsSum += edge->weight;
    graph->matrix[edge->n1][edge->n2] = edge->weight;
    graph->matrix[edge->n2][edge->n1] = edge->weight;
}

/**
 * @brief Removed edge from graph
 * 
 * @param graph given graph
 * @param edge edge to be removed
 */
void removeEdge(GraphStruct* graph, Edge* edge) {
    graph->edgesNum -= 1;
    graph->weightsSum -= edge->weight;
    graph->matrix[edge->n1][edge->n2] = 0;
    graph->matrix[edge->n2][edge->n1] = 0;
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
            if (graph->matrix[i][j] != 0) {
                edges[edgeID] = new Edge;
                edges[edgeID]->isUsed = 0;
                edges[edgeID]->weight = graph->matrix[i][j];
                edges[edgeID]->n1 = i;
                edges[edgeID]->n2 = j;
                edgeID++;
            }
        }
    }

    
    // sort the edges - buble sort descending
    Edge * tmpEdge = NULL;
    for (int i = 0; i < graph->edgesNum; i++) {
        tmpEdge = edges[i];
        for (int j = i+1; j < graph->edgesNum; j++) {
            if (edges[j]->weight > edges[i]->weight) {
                edges[i] = edges[j];
                edges[j] = tmpEdge;
                // restart the highest
                tmpEdge = edges[i];
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
void addResult(Results* results, GraphStruct* graph, int* cNodes) {
    // delete old worse solutions
    if (results->results != NULL) {
        if (results->results->graph->weightsSum < graph->weightsSum) {
            ResultNode* tmpResult = results->results;
            while (tmpResult->next != NULL) {
                tmpResult = results->results->next;
                delete results->results;
                results->results = NULL;
                results->results = tmpResult;
            }
            delete tmpResult;
            results->results = NULL;
        }
    }

    // add new result
    ResultNode* newResult = new ResultNode;
    newResult->cNodes = copyVector(cNodes, graph->nodesNum);
    newResult->graph = copyGraph(graph);
    newResult->next = results->results;

    results->results = newResult;
}