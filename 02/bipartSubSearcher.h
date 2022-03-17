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

#include "graph.h"

void searchColorCombination(Graph* graph, Graph* subgraph, Results* results, int* cNodes, int trashWeights, int lastEdgeId, Edge* edge, int cN1, int cN2);

/**
 * @brief Function combinatons of edges 0-1, 1-0 or none
 * 
 * @param graph given graph to be used
 * @param subgraph subgraph with result
 * @param results all best results
 * @param cNodes nodes colors
 * @param trashWeights weight of edges that cannot be used
 */
void searchBiCoSubgraphs(Graph* graph, Graph* subgraph, Results* results, int* cNodes, int trashWeights, int lastEdgeId) {
    // test if should continue
    if (results->results != NULL) {
        int freeWeights = graph->getWeightsSum() - trashWeights - subgraph->getWeightsSum();
       
        if ( (subgraph->getWeightsSum() + freeWeights) < results->results->weight ) {
            delete subgraph;
            delete[] cNodes;
            return;
        };
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
        if (results->results != NULL && results->results->weight > subgraph->getWeightsSum()) {
            delete subgraph;
            delete[] cNodes;
            return;
        };

        // subgraph is not valid
        if (isConnected(subgraph) == false) {
            delete subgraph;
            delete[] cNodes;
            return;
        };
     
        // store results
        #pragma omp critical
            addResult(results, subgraph, cNodes, graph->getEdgesNum());
        return;
    }
    

    // if can be colored edge 0 1 and delete edge after you are done
    if (
        cNodes[edge->n1] != 1 && cNodes[edge->n2] != 0 &&
        canColorNode(subgraph, cNodes, edge->n1, 0) &&
        canColorNode(subgraph, cNodes, edge->n2, 1)
        ) {
            #pragma omp task
                searchColorCombination(graph, subgraph->copyGraph(), results, copyVector(cNodes, graph->getNodesNum()), trashWeights, curEdgeId+1, edge, 0, 1);
            //#pragma omp taskwait
    }

    // if can be colored edge 1 0 and delete edge after you are done
    if (
        cNodes[edge->n1] != 0 && cNodes[edge->n2] != 1 &&
        canColorNode(subgraph, cNodes, edge->n1, 1) &&
        canColorNode(subgraph, cNodes, edge->n2, 0)
        ) {
            #pragma omp task
                searchColorCombination(graph, subgraph->copyGraph(), results, copyVector(cNodes, graph->getNodesNum()), trashWeights, curEdgeId+1, edge, 1, 0);
            //#pragma omp taskwait
    }

    // dont use this edge
    trashWeights += edge->weight;
    #pragma omp task
        searchBiCoSubgraphs(graph, subgraph->copyGraph(), results, copyVector(cNodes, graph->getNodesNum()), trashWeights, curEdgeId+1);
    #pragma omp taskwait
    delete subgraph;
    delete[] cNodes;
}

/**
 * @brief todo
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
void searchColorCombination(Graph* graph, Graph* subgraph, Results* results, int* cNodes, int trashWeights, int lastEdgeId, Edge* edge, int cN1, int cN2) {
        // color the nodes
        cNodes[edge->n1] = cN1;
        cNodes[edge->n2] = cN2;

        // add the edge to the graph
        subgraph->addEdge(edge);

        // continue the search with this setting
        searchBiCoSubgraphs(graph, subgraph, results, cNodes, trashWeights, lastEdgeId);
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

    // parallel run
    #pragma omp parallel 
    {
        // get the results
        #pragma omp single
            searchBiCoSubgraphs(graph, subgraph, results, cNodes, 0, 0);
    }

    // delete unnecesary
    //delete[] cNodes;
    //for (int i = 0; i < graph->getEdgesNum(); i++) delete graph->getEdges()[i];
    //delete[] graph->getEdges();
    //delete subgraph;

    // return the result
    return results;
}
