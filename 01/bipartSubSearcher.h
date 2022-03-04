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

#include "graph.h"

/**
 * @brief Function combinatons of edges 0-1, 1-0 or none
 * 
 * @param graph given graph to be used
 * @param subgraph subgraph with result
 * @param edges set of all edges
 * @param results all best results
 * @param cNodes nodes colors
 * @param trashWeights weight of edges that cannot be used
 */
void searchBiCoSubgraphs(GraphStruct* graph, GraphStruct* subgraph, Edge** edges, Results* results, int* cNodes, int trashWeights) {
    // test if should continue
    if (results->results != NULL) {
        int freeWeights = graph->weightsSum - trashWeights - subgraph->weightsSum;
       
        if ( (subgraph->weightsSum + freeWeights) < results->results->graph->weightsSum ) return;
    }

    // get first available edge
    Edge* edge = NULL;
    for (int i = 0; i < graph->edgesNum; i++) {
        if (edges[i]->isUsed == 0) {
            edge = edges[i];
            break;
        }
    }

    // no more available edge
    if (edge == NULL) {
        // Is new solution better
        if (results->results != NULL && results->results->graph->weightsSum > subgraph->weightsSum) return;

        // subgraph is not valid
        if (isConnected(subgraph) == false) return;
     
        // store results
        addResult(results, subgraph, cNodes);
        return;
    }
    

    // if can be colored edge 0 1 and delete edge after you are done
    if (
        cNodes[edge->n1] != 1 && cNodes[edge->n2] != 0 &&
        canColorNode(subgraph, cNodes, edge->n1, 0) &&
        canColorNode(subgraph, cNodes, edge->n2, 1)
        ) {
        // edge is being used
        edge->isUsed = 1;
        
        // indication to clean after search
        int n1C = cNodes[edge->n1] == 0 ? 0 : 1;
        int n2C = cNodes[edge->n2] == 1 ? 0 : 1;

        // color the nodes
        cNodes[edge->n1] = 0;
        cNodes[edge->n2] = 1;

        // add the edge to the graph
        addEdge(subgraph, edge);

        // continue the search with this setting
        searchBiCoSubgraphs(graph, subgraph, edges, results, cNodes, trashWeights);

        // remove edge from graph
        removeEdge(subgraph, edge);

        cNodes[edge->n1] = n1C == 1 ? -1 : 0;
        cNodes[edge->n2] = n2C == 1 ? -1 : 1;
    }

    // if can be colored edge 1 0 and delete edge after you are done
    if (
        cNodes[edge->n1] != 0 && cNodes[edge->n2] != 1 &&
        canColorNode(subgraph, cNodes, edge->n1, 1) &&
        canColorNode(subgraph, cNodes, edge->n2, 0)
        ) {
        // edge is being used
        edge->isUsed = 1;
        
        // indication to clean after search
        int n1C = cNodes[edge->n1] == 1 ? 0 : 1;
        int n2C = cNodes[edge->n2] == 0 ? 0 : 1;

        // color the nodes
        cNodes[edge->n1] = 1;
        cNodes[edge->n2] = 0;

        // add the edge to the graph
        addEdge(subgraph, edge);

        // continue the search with this setting
        searchBiCoSubgraphs(graph, subgraph, edges, results, cNodes, trashWeights);

        // remove edge from graph
        removeEdge(subgraph, edge);

        cNodes[edge->n1] = n1C == 1 ? -1 : 1;
        cNodes[edge->n2] = n2C == 1 ? -1 : 0;
    }

    // dont use this edge
    edge->isUsed = -1;
    trashWeights += edge->weight;
    searchBiCoSubgraphs(graph, subgraph, edges, results, cNodes, trashWeights);
    
    edge->isUsed = 0;
}

Results* getMaxBiparSubgraph(std::string graphName) {
    GraphStruct* graph = loadGraph(graphName);

    if (graph == NULL) return NULL;

    // test if graph can be used
    if (!isConnected(graph)) {
        std::cout << "Graph is not connected" << '\n';
        return 0;
    }
  
    // end if the graph itself is biparted
    if (isBiparted(graph)) {
        std::cout << "Graph is biparted" << '\n';
    }

    // create all neccesary components
    Results* results = new Results;
    results->results = NULL;

    GraphStruct* subgraph = createGraph(graph->nodesNum);
    Edge** edges = createSorEdgesLL(graph);
    
    int* cNodes = new int[graph->nodesNum];
    for (int i = 0; i < graph->nodesNum; i++) {
        cNodes[i] = -1;
    }

    // color one node
    cNodes[edges[0]->n1] = 0;

    // get the results
    searchBiCoSubgraphs(graph, subgraph, edges, results, cNodes, 0);

    // delete unnecesary
    delete[] cNodes;
    for (int i = 0; i < graph->edgesNum; i++) delete edges[i];
    delete[] edges;
    delete subgraph;
    delete graph;

    // return the result
    return results;
}
