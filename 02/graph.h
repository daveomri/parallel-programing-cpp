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
#include <vector>
#include <numeric>

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

    Edge(int weight, int n1, int n2, int isUsed): weight(weight), n1(n1), n2(n2), isUsed(isUsed) {}

    Edge(Edge* edge) {
        this->weight = edge->weight;
        this->n1 = edge->n1; 
        this-> n2 = edge-> n2; 
        this->isUsed = edge->isUsed; 
    }
};

/**
 * @brief Class represents the graph
 * 
 */
class Graph {
    private:
        int nodesNum;
        int weightsSum;
        int edgesNum;
        int* matrix;
    public:
    //--------------------------------------
        Graph(int numNodes): nodesNum(numNodes), weightsSum(0), edgesNum(0) {
            matrix = new int[numNodes*numNodes];
            for (int i = 0; i < numNodes*numNodes; i++) {
                matrix[i] = 0;
            }
        }

        Graph(std::string graphName) {
            this->nodesNum = 0;
            this->weightsSum = 0;
            this->edgesNum = 0;
            loadGraph(graphName);
        }

        ~Graph() {
            delete[] matrix;
        }

        int getNodesNum() {return this->nodesNum;}

        int getWeightsSum() {return this->weightsSum;}

        int getEdgesNum() {return this->edgesNum;}

        void setWeightsSum(int weightsSum) {
            this->weightsSum = weightsSum;
        }

        void setEdgesNum(int edgesNum) {
            this->edgesNum = edgesNum;
        }

        void setEdgeWeight(int x, int y, int weight) {
            if (matrix == NULL) return;
            matrix[x*this->nodesNum + y] = weight;
        }

        int getEdgeWeight(int x, int y) {
            if (matrix == NULL) return -1;
            return matrix[x*this->nodesNum + y];
        }

        void copyGraph(Graph* graph) {
            int matrixSize = graph->getNodesNum()*graph->getNodesNum();
            for (int i = 0; i < matrixSize; i++) {
                this->matrix[i] = graph->matrix[i];
            }
            this->nodesNum = graph->getNodesNum();
            this->edgesNum = graph->getEdgesNum();
            this->weightsSum = graph->weightsSum;
        }

        /**
         * @brief Function copies given graph
         * 
         * @param graph given graph
         * @return Graph* new identical graph
         */
        Graph* copyGraph() {
            Graph* newGraph = new Graph(this->getNodesNum());
            
            newGraph->copyGraph(this);
            
            return newGraph;
        }

        /**
         * @brief Method loads the file with given name
         *        in the graph_mbp folder and stores its
         *        data to the array
         * 
         * @param fileName File name with graph prescription
         */
        void loadGraph(std::string graphName) {
            // load the file
            std::string line;
            std::ifstream graphFile ("graf_mbp/" + graphName);
            
            if (!graphFile.is_open()) {
                std::cout << "Unable to open graf_mbp/" + graphName + "\n";
                return;
            }

            // read the first line
            int num;
            if (getline (graphFile, line) ) {
                std::istringstream iss(line);

                iss >> num;
            }
            else std::cout << "Unable to read the first line.\n";
            
            // stop method if input matrix dimension is bad
            if (num <= 0) return;

            // create matrix
            this->matrix = new int[num*num];
            this->nodesNum = num;
            
            int tmpWeight = 0;

            // read the rest of the file
            for (int i = 0; i < num; i++) {
                std::getline(graphFile, line);
                std::istringstream iss(line);

                for (int j = 0; j < num; j++) {
                    iss >> tmpWeight;
                    this->setEdgeWeight(i, j, tmpWeight);
                    if (tmpWeight != 0){
                        this->edgesNum++;
                        this->weightsSum+=tmpWeight;
                    }
                }
            }

            // fix the total weight and total edges num
            this->weightsSum = this->weightsSum/2;
            this->edgesNum = this->edgesNum/2;
        }

        /**
         * @brief Function adds edge to graph
         * 
         * @param edge edge to be added
         */
        void addEdge(Edge* edge) {

            this->edgesNum += 1;
            this->weightsSum += edge->weight;
            this->setEdgeWeight(edge->n1, edge->n2, edge->weight);
            this->setEdgeWeight(edge->n2, edge->n1, edge->weight);
        }

        /**
         * @brief Removed edge from graph
         * 
         * @param edge edge to be removed
         */
        void removeEdge(Edge* edge) {
            this->edgesNum-=1;
            this->weightsSum -= edge->weight;
            this->setEdgeWeight(edge->n1, edge->n2, 0);
            this->setEdgeWeight(edge->n2, edge->n1, 0);
        }

        /**
         * @brief Function sums the edges weights
         * 
         * @return int 
         */
        int sumEdgesWeights(){
            int edgeSum = 0;
            for (int i = 0; i < this->getNodesNum(); i++) {
                for (int j = i; j < this->getNodesNum(); j++) {
                    edgeSum += this->getEdgeWeight(i, j);
                }
            }
            return edgeSum;
        }
};

// ----------------------------------------------

// Additional graph functions

/**
 * @brief Structure represents the search result
 * 
 */
struct ResultNode {
    ResultNode* next;

    // data
    //Graph* graph;
    Edge** edges;
    int edgesNum;
    int* cNodes;
    int weight;
    ~ResultNode() {
        for (int i = 0; i < edgesNum; i++) {
            delete edges[i];
        }
        delete[] edges;
        delete [] cNodes;
    }
};

struct Results
{
    ResultNode* results;
};

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
bool canBeColored(Graph* graph, int* cNodes, int cNode, int color) {
    for (int i = 0; i < graph->getNodesNum(); i++) {
        if (graph->getEdgeWeight(cNode, i) != 0 && cNodes[i] == color) return false;
    }
    return true;
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
bool canColorNode(Graph* graph, int* cNodes, int cNode, int color) {
    for (int i = 0; i < graph->getNodesNum(); i++) {
        if (graph->getEdgeWeight(cNode, i) != 0 && cNodes[i] == color) return false;
    }
    return true;
}

/**
 * @brief Create array of sorted edges
 * 
 * @param graph graph with edges to be sorted
 * @return Edge** sorted array of edges
 */
Edge** createSorEdgesLL(Graph* graph) {
    Edge** edges = new Edge*[graph->getEdgesNum()];
    int edgeID = 0;

    // add the edges
    for (int i = 0; i < graph->getNodesNum(); i++) {
        for (int j = i; j < graph->getNodesNum(); j++) {
            if (graph->getEdgeWeight(i, j) != 0) {
                edges[edgeID] = new Edge(graph->getEdgeWeight(i,j), i, j, 0);
                edgeID++;
            }
        }
    }

    
    // sort the edges - buble sort descending
    Edge * tmpEdge = NULL;
    for (int i = 0; i < graph->getEdgesNum(); i++) {
        tmpEdge = edges[i];
        for (int j = i+1; j < graph->getEdgesNum(); j++) {
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
void srchdfs(Graph* graph, int* visNodes, int curNode) {
    visNodes[curNode] = 1;
    for (int i = 0; i < graph->getNodesNum(); i++) {
        if ( visNodes[i] == 0 && graph->getEdgeWeight(curNode, i) != 0) {
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
bool isConnected(Graph* graph) {
    // create graph unvisited nodes
    int* visNodes = new int[graph->getNodesNum()];
    
    // sets zeros to all nodes
    for (int i = 0; i < graph->getNodesNum(); i++) {
        visNodes[i] = 0;
    } 

    // search graph - start with 0
    srchdfs(graph, visNodes, 0);

    // test the visited nodes
    for (int i = 0; i < graph->getNodesNum(); i++) {
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
bool colorGraph(Graph* graph, int* cNodes, int cNode, int color) {
    // check if the color can be used
    if (canBeColored(graph, cNodes, cNode, color) == false) return false;

    // color the node
    cNodes[cNode] = color;

    int nextColor = color == 0 ? 1 : 0;

    // continue the coloring
    for ( int i = 0; i  < graph->getNodesNum(); i++ ) {
        if (graph->getEdgeWeight(cNode, i) != 0 && cNodes[i] == -1) {
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
bool isBiparted(Graph* graph) {
    // create graph unvisited nodes
    int* visNodes = new int[graph->getNodesNum()];
    for (int i = 0; i < graph->getNodesNum(); i++) {
        // set the color to undefined
        visNodes[i] = -1;
    }

    if (colorGraph(graph, visNodes, 0, 0)) {
        for (int i = 0; i < graph->getNodesNum(); i++) {
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

Edge** copyEdges(Edge** edges, int sizeEdges, int sizeNew) {
    Edge** newEdges = new Edge*[sizeNew]; 
    int tmpId = 0;

    for (int i = 0; i < sizeEdges; i++) {
        if (edges[i]->isUsed == 1) {
            newEdges[tmpId++] = new Edge(edges[i]);
        }
    }
    return newEdges;
}

/**
 * @brief Function adds new result to top od the linked list
 * 
 * @param results results so far, linked list
 * @param graph graph with result
 * @param cNodes colored nodes of given graph1
 * @return ResultNode* updated results
 */
void addResult(Results* results, Graph* graph, int* cNodes, Edge** edges, int edgesNum) {
    // delete old worse solutions
    if (results->results != NULL) {
        if (results->results->weight < graph->getWeightsSum()) {
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
    newResult->cNodes = copyVector(cNodes, graph->getNodesNum());
    //newResult->graph = graph->copyGraph();
    newResult->edges = copyEdges(edges, edgesNum, graph->getEdgesNum());
    newResult->edgesNum = graph->getEdgesNum();
    newResult->weight = graph->getWeightsSum();
    newResult->next = results->results;

    results->results = newResult;
}