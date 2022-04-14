/**
 * @file bipartSubSearcher.h
 * @author David Omrai
 * @brief All necessary functions for bipart. max subgraph search
 * @version 0.1
 * @date 2022-03-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <omp.h>
#include <list>
#include <mpi.h>

#include "graph.h"
#define THREADS 6

/**
 * @brief Structure represents state of search
 * 
 */
struct SearchState {
    Graph* subgraph=NULL;
    int* cNodes=NULL;
    int trashWeights=0;
    int lastEdgeId=0;
    //-------------
    SearchState(Graph*subgraph,int*cNodes,int trashWeights,int lastEdgeId):
        subgraph(subgraph),cNodes(cNodes),trashWeights(trashWeights),lastEdgeId(lastEdgeId) {

    }
    ~SearchState() {
        if (this->subgraph) delete subgraph;
        //if (this->cNodes) delete[] cNodes;
    }
};

struct Message {
    int trashWeights;
    int lastEdgeId;
    int isEnd;
    int procId;
};

void setColorCombination(SearchState* searchState, Edge* edge, int cN1, int cN2);


/**
 * @brief Funciton uses BFS approach to generate given number of search states
 * 
 * @param graph 
 * @param results 
 * @param list 
 * @param statesNum 
 */
void bfsSearchStates(Graph* graph, Results* results, std::list<SearchState*> &list, size_t statesNum) {
    while(list.size() < statesNum) {
        // End search if nothing remained
        if (list.size() == 0) return;

        // Get new state
        SearchState* searchState = list.back();
        list.pop_back();
        
        // Prepare variables
        Graph* subgraph = searchState->subgraph;
        int* cNodes = searchState->cNodes;
        int trashWeights = searchState->trashWeights;
        int lastEdgeId = searchState->lastEdgeId;
        
        // Test if to continue
        if (results->results != NULL) {
            int freeWeights = graph->getWeightsSum() - trashWeights - subgraph->getWeightsSum();
        
            if ( (subgraph->getWeightsSum() + freeWeights) < results->results->weight ) {
                delete searchState;
                continue;
            }
        }

        // get first available edge
        Edge* edge = NULL;
        int curEdgeId = 0;
        //int a;
        for (int i = lastEdgeId; i < graph->getEdgesNum(); i++) {
            if (!subgraph->getEdgeWeight(graph->getEdges()[i]->n1, graph->getEdges()[i]->n2)) {
                edge = graph->getEdges()[i];
                curEdgeId = i;
                break;
            }
        }

        // no more available edge
        if (edge == NULL) {
            
            // Is new solution better
            if ((results->results == NULL && isConnected(subgraph) == true) || 
                (results->results != NULL && results->results->weight <= subgraph->getWeightsSum() && isConnected(subgraph) == true)) {
                addResult(results, subgraph, cNodes);
            }
            
            delete searchState;
            continue;
        }

        // if can be colored edge 0 1 and delete edge after you are done
        if (
            cNodes[edge->n1] != 1 && cNodes[edge->n2] != 0 &&
            canColorNode(subgraph, cNodes, edge->n1, 0) &&
            canColorNode(subgraph, cNodes, edge->n2, 1)
            ) {
                Graph* newSubgraph = subgraph->copyGraph();
                int* newCNodes = copyVector(cNodes, graph->getNodesNum());
                SearchState* newState = new SearchState(newSubgraph, newCNodes, trashWeights, curEdgeId+1);
                setColorCombination(newState, edge, 0, 1);
                // continue the search with this setting
                list.push_front(newState);
        }

        // if can be colored edge 1 0 and delete edge after you are done
        if (
            cNodes[edge->n1] != 0 && cNodes[edge->n2] != 1 &&
            canColorNode(subgraph, cNodes, edge->n1, 1) &&
            canColorNode(subgraph, cNodes, edge->n2, 0)
            ) {
                Graph* newSubgraph = subgraph->copyGraph();
                int* newCNodes = copyVector(cNodes, graph->getNodesNum());
                SearchState* newState = new SearchState(newSubgraph, newCNodes, trashWeights, curEdgeId+1);
                setColorCombination(newState, edge, 1, 0);
                // continue the search with this setting
                list.push_front(newState);
        }

        // dont use this edge
        trashWeights += edge->weight;
        Graph* newSubgraph = subgraph->copyGraph();
        int* newCNodes = copyVector(cNodes, graph->getNodesNum());
        list.push_front(new SearchState(newSubgraph, newCNodes, trashWeights, curEdgeId+1));
        // clean the mess
        delete searchState;
    }
}

/**
 * @brief Function combinatons of edges 0-1, 1-0 or none
 * 
 * @param graph given graph to be used
 * @param subgraph subgraph with result
 * @param results all best results
 * @param cNodes nodes colors
 * @param trashWeights weight of edges that cannot be used
 */
void searchBiCoSubgraphs(Graph* graph, Results* results, SearchState* searchState) {
    // Define search data
    Graph* subgraph = searchState->subgraph;
    int* cNodes = searchState->cNodes;
    int trashWeights = searchState->trashWeights;
    int lastEdgeId = searchState->lastEdgeId;


    // test if should continue
    if (results->results != NULL) {
        int freeWeights = graph->getWeightsSum() - trashWeights - subgraph->getWeightsSum();
       
        if ( (subgraph->getWeightsSum() + freeWeights) < results->results->weight ) {
            delete searchState;
            return;
        }
    }

    // get first available edge
    Edge* edge = NULL;
    int curEdgeId = 0;
    //int a;
    for (int i = lastEdgeId; i < graph->getEdgesNum(); i++) {
        if (!subgraph->getEdgeWeight(graph->getEdges()[i]->n1, graph->getEdges()[i]->n2)) {
            edge = graph->getEdges()[i];
            curEdgeId = i;
            break;
        }
    }

    // no more available edge
    if (edge == NULL) {
        // Is new solution better
        if ((results->results == NULL && isConnected(subgraph) == true) || 
            (results->results != NULL && results->results->weight <= subgraph->getWeightsSum() && isConnected(subgraph) == true)) {
            addResult(results, subgraph, cNodes);
        
        }
        delete searchState;
        return;
    }
    

    // if can be colored edge 0 1 and delete edge after you are done
    if (
        cNodes[edge->n1] != 1 && cNodes[edge->n2] != 0 &&
        canColorNode(subgraph, cNodes, edge->n1, 0) &&
        canColorNode(subgraph, cNodes, edge->n2, 1)
        ) {
            Graph* newSubgraph = subgraph->copyGraph();
            int* newCNodes = copyVector(cNodes, graph->getNodesNum());
            SearchState* newState = new SearchState(newSubgraph, newCNodes, trashWeights, curEdgeId+1);
            setColorCombination(newState, edge, 0, 1);
            // continue the search with this setting
            searchBiCoSubgraphs(graph, results, newState);
    }

    // if can be colored edge 1 0 and delete edge after you are done
    if (
        cNodes[edge->n1] != 0 && cNodes[edge->n2] != 1 &&
        canColorNode(subgraph, cNodes, edge->n1, 1) &&
        canColorNode(subgraph, cNodes, edge->n2, 0)
        ) {
            Graph* newSubgraph = subgraph->copyGraph();
            int* newCNodes = copyVector(cNodes, graph->getNodesNum());
            SearchState* newState = new SearchState(newSubgraph, newCNodes, trashWeights, curEdgeId+1);
            setColorCombination(newState, edge, 1, 0);
            // continue the search with this setting
            searchBiCoSubgraphs(graph, results, newState);
    }

    // dont use this edge
    trashWeights += edge->weight;
    Graph* newSubgraph = subgraph->copyGraph();
    int* newCNodes = copyVector(cNodes, graph->getNodesNum());
    searchBiCoSubgraphs(graph, results, new SearchState(newSubgraph, newCNodes, trashWeights, curEdgeId+1));
    // clean the mess
    delete searchState;
}


void copyArray(int** narr, int** parr, int size) {
    for (int i = 0; i < size; i++) {
        narr[i] = parr[i];
    }
}

Graph* createSubgraph(Graph* graph, int* cNodes) {
    // create subgraph
    Graph* subgraph = new Graph(graph->getNodesNum());

    for (int i = 0; i < graph->getNodesNum(); i++) {
        for (int j = i+1; j < graph->getNodesNum(); j++) {
            if (cNodes[i] != cNodes[j] && cNodes[i] != -1 && cNodes[j] != -1 && graph->getEdgeWeight(i, j) != 0) {
                Edge* newEdge = new Edge(graph->getEdgeWeight(i, j), i, j, 1);
                subgraph->addEdge(newEdge);
                delete newEdge;
            }
        }
    }

    return subgraph;
}

void searchBiCoSubgraphsParallel(int &argc, char **argv) {
    int numProcs;
    int rank;
    char name[80];
    int length;
    const int tag = 13;

    // Initialize the MPI
    MPI_Init(&argc, &argv);
    // Number of ranks in each communicator
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    // Which proces am I - mpi_comm_world is comunicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the name of the processor
    MPI_Get_processor_name(name, &length);

    // struct creation for state
    const int nitems = 4;
    int blockLens[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype mpi_state_type;
    MPI_Aint offsets[4];

    offsets[0] = offsetof(Message, trashWeights);
    offsets[1] = offsetof(Message, lastEdgeId);
    offsets[2] = offsetof(Message, isEnd);
    offsets[3] = offsetof(Message, procId);

    MPI_Type_create_struct(nitems, blockLens, offsets, types, &mpi_state_type);
    MPI_Type_commit(&mpi_state_type); 


    // The main part ----------------------------------------------------------
    if (rank == 0) {
        // ---------------------------------------
        std::string graphName = argv[1];
        Graph* graph = new Graph(graphName);
        // ---------------------------------------

        // create all neccesary components
        Results* results = new Results;
        results->results = NULL;

        Graph* subgraph = new Graph(graph->getNodesNum());
        graph->setEdges(createSorEdgesLL(graph));

        
        int* cNodes = new int[graph->getNodesNum()];
        for (int i = 0; i < graph->getNodesNum(); i++) {
            cNodes[i] = -1;
        }

        // color one node
        cNodes[graph->getEdges()[0]->n1] = 0;

        // BFS search to generate enough search states
        std::list<SearchState*> searchStates;
        searchStates.push_back(new SearchState(subgraph, cNodes, 0, 0));
        bfsSearchStates(graph, results, searchStates, 3);

        // ---------------------------------------

        // cycle through the list
        MPI_Status status;
        Message recv;

        int curSlave = 0;
        while (!searchStates.empty()) {
            printf("list size is %d \n", (int)searchStates.size());
            // choose which proces to work with
            curSlave+=1;
            if (curSlave>= 2) {
                curSlave = 1;
            }
            printf("Current slave is %d\n", curSlave);
            // --------------------------------
            // receive request for work
            MPI_Recv(&recv, 1, mpi_state_type, curSlave, tag, MPI_COMM_WORLD, &status);
            printf("just received from %d\n", recv.procId);
            // give job --------------------------------------
            Message send;

            SearchState* srchsSate = searchStates.front();
            searchStates.pop_front();
            // ---------------------------------------
            send.lastEdgeId = srchsSate->lastEdgeId;
            send.trashWeights = srchsSate->trashWeights;
            //send.cNodes = srchsSate->cNodes;
            send.isEnd = 0;
            send.procId = recv.procId;
            // ---------------------------------------

            const int dest = recv.procId;
            MPI_Send(&send, 1, mpi_state_type, dest, 6, MPI_COMM_WORLD);
            
            printf("this is cnodes %d\n", srchsSate->cNodes[0]);
            // MPI_Send(srchsSate->cNodes, graph->getNodesNum(), MPI_INT, dest, 9, MPI_COMM_WORLD);
            // delete srchsSate;
            // ---------------------------------------------------------------
        }
        printf("endos\n\n\n\n");
        // end all proceses
        for (int i = 1; i < 2; i++) {
            // wait till comes the response
            MPI_Recv(&recv, 1, mpi_state_type, i, tag, MPI_COMM_WORLD, &status);

            const int dest = recv.procId;
            Message send;
            // ---------------------------------------
            send.lastEdgeId = 0;
            send.trashWeights = 0;
            send.isEnd = 1;
            send.procId = recv.procId;
            // ---------------------------------------
            MPI_Send(&send, 1, mpi_state_type, dest, 6, MPI_COMM_WORLD);
        }
        printf("done first half\n");
        // get results from processes


        for (int i = 1; i < 2; i++) {
            int ccNodes;
            printf("this far");
            MPI_Recv(&ccNodes, 1, MPI_INT, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //Graph* subgraph = createSubgraph(graph, cNodes);
            // todo multiple results - while messages keeps on comming
            //addResult(results, subgraph, cNodes);
            //delete subgraph;
        }

        // print results here
        printf("here and done\n");
        if (false) { //results->results != NULL
            std::cout << "RESULTS" << std::endl;
            
            ResultNode* tmpRes = results->results;
            while (tmpRes != NULL) {
                std::cout << "Weight sum: " << std::endl;
                std::cout << tmpRes->weight << std::endl;
                std::cout << "Node colors:" << std::endl;
                for (int i = 0; i < graph->getNodesNum(); i++) {
                    std::cout << i << ":" << tmpRes->cNodes[i] << ", ";
                }
                std::cout << '\n';

                // cout << "Edges:" << endl;
                // for (int i = 0; i < tmpRes->graph->getNodesNum(); i++) {
                //     for (int j = 0; j < tmpRes->graph->getNodesNum(); j++) {
                //         cout << " | " << setw(4) << tmpRes->graph->getEdgeWeight(i, j);
                //     }
                //     cout << " |" << endl;
                // }

                std::cout << "Edges:" << "\n";
                int mi = 0;
                for (int i = 0; i < graph->getNodesNum(); i++) {
                    for (int j = 0; j < graph->getNodesNum(); j++) {
                        mi = graph->getMatrixPosition(i, j);
                        if (mi == -1) std::cout << "  |   " << std::setw(4) << 0;
                        else std::cout << "  |   " << std::setw(4) << tmpRes->graphWeights[mi];
                        
                    }
                    std::cout << "  |" << "\n";
                }

                tmpRes = tmpRes->next;
                delete results->results;
                results->results = tmpRes;

                std::cout << "-----------------------------------------------------------" << '\n';
            }
        }
        delete results;
        delete graph;
    }
    else if (rank == 1) {
        MPI_Status status;
        // Create for each process-----------
        std::string graphName = argv[1];
        Graph* graph = new Graph(graphName);
        // Get results
        Results* results = new Results;
        results->results = NULL;
        //-----------------------------------

        bool END = false; 
        // Cycle-----------------------------
        while (!END) {
            // Send to get some work
            Message send;
            send.procId = rank;
            printf("Imma waiting %d\n", rank);
            MPI_Send(&send, 1, mpi_state_type, 0, tag, MPI_COMM_WORLD);
            printf("Imma working %d\n", rank);

            // Wait for the response
            Message recv;
            
            MPI_Recv(&recv, 1, mpi_state_type, 0, 6, MPI_COMM_WORLD, &status);
            
            
            // check if it is the end ------------------------------------------
            if (recv.isEnd) {
                printf("------------This is the end\n");
                END = true;
            }
            else {
                // int* ccNodes;
                // MPI_Recv(ccNodes, graph->getNodesNum(), MPI_INT, 0, 9, MPI_COMM_WORLD, &status);
                // printf("this is cnodes slave %d\n", ccNodes[0]);
                // printf("here\n");

                // printf("just received slave data %d %d %d %d\n", recv.isEnd, recv.lastEdgeId, recv.procId, recv.trashWeights);
                // // Create the graph --------------------------------------------
                // printf("There could be error\n");
                // Graph* subgraph = createSubgraph(graph, ccNodes);
                // printf("yup there it is\n");
                // SearchState* searchState = new SearchState(subgraph, ccNodes, recv.trashWeights, recv.lastEdgeId);
                // printf("Nope it's here\n");
                // // Find best solution --- possibly paralel
                // //searchBiCoSubgraphs(graph, results, searchState);
                // delete subgraph;
            }
        }

        // return the result
       
        // int* cNodes;
        // if (results->results) {
        //     printf("1\n");
        //     cNodes = results->results->cNodes;
        // }
        // else {
        //     printf("2\n");
        //     cNodes = new int[graph->getNodesNum()];
        //     printf("here\n");
        //     cNodes[0] = -1;
        // }
        printf("yet working\n");
        int cccNodes = 99;
        MPI_Send(&cccNodes, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        printf("this is bad or what\n");
        delete results;
        delete graph;
    }

    MPI_Type_free(&mpi_state_type);
    MPI_Finalize();
    printf("here am i \n\n");
}

/**
 * @brief Function returns sorted array of states by their total weight
 * 
 * @param list 
 * @return SearchState** 
 */
SearchState** sortList(std::list<SearchState*>&list) {
    SearchState**sorArr = new SearchState*[list.size()];
    int lastIndex = 0;
    for (auto state: list) {
        sorArr[lastIndex++] = state;
    }
    return sorArr;
}

/**
 * @brief Function sets the color of nodes and then continue the search
 * 
 * @param graph 
 * @param subgraph 
 * @param results 
 * @param cNodes 
 * @param trashWeights 
 * @param lastEdgeId 
 * @param cN1 
 * @param cN2 
 */
void setColorCombination(SearchState* searchState, Edge* edge, int cN1, int cN2) {
        // color the nodes
        searchState->cNodes[edge->n1] = cN1;
        searchState->cNodes[edge->n2] = cN2;

        // add the edge to the graph
        searchState->subgraph->addEdge(edge);
}

void getMaxBiparSubgraph(int &argc, char **argv) {

    // if (graph == NULL) return;

    // // test if graph can be used
    // if (!isConnected(graph)) {
    //     std::cout << "Graph is not connected" << '\n';
    //     return;
    // }

    // // end if the graph itself is biparted
    // if (isBiparted(graph)) {
    //     std::cout << "Graph is biparted" << '\n';
    //     return;
    // }

    searchBiCoSubgraphsParallel(argc, argv);
}
