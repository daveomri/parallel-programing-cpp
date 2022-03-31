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
        if (this->cNodes) delete[] cNodes;
    }
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
                addResult(results, subgraph, cNodes, graph->getEdgesNum());
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
        #pragma omp critical 
        {
            if ((results->results == NULL && isConnected(subgraph) == true) || 
                (results->results != NULL && results->results->weight <= subgraph->getWeightsSum() && isConnected(subgraph) == true)) {
                addResult(results, subgraph, cNodes, graph->getEdgesNum());
            
            }
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

Results* getMaxBiparSubgraph(Graph* graph) {

    if (graph == NULL) return NULL;

    // test if graph can be used
    if (!isConnected(graph)) {
        std::cout << "Graph is not connected" << '\n';
        return NULL;
    }
  
    // end if the graph itself is biparted
    if (isBiparted(graph)) {
        std::cout << "Graph is biparted" << '\n';
    }

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
    std::list<SearchState*> list;
    list.push_back(new SearchState(subgraph, cNodes, 0, 0));
    bfsSearchStates(graph, results, list, 400);
    SearchState**states = sortList(list);

    // for cycle for paralel cycle
    omp_set_num_threads(THREADS);
    int i = 0;
   
    #pragma omp parallel for schedule(guided) num_threads(THREADS)
    for (i = 0; i < (int)list.size(); i++ ) {
        //std::cout << i << "\n";
        searchBiCoSubgraphs(graph, results, states[i]);
    }
    
    // Clean the mess
    delete states;

    // return the result
    return results;
}
