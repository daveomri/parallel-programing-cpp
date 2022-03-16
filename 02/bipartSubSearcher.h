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
       
        if ( (subgraph->getWeightsSum() + freeWeights) < results->results->weight ) return;
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
        if (results->results != NULL && results->results->weight > subgraph->getWeightsSum()) return;

        // subgraph is not valid
        if (isConnected(subgraph) == false) return;
     
        // store results
        addResult(results, subgraph, cNodes, graph->getEdgesNum());
        return;
    }
    

    // if can be colored edge 0 1 and delete edge after you are done
    if (
        cNodes[edge->n1] != 1 && cNodes[edge->n2] != 0 &&
        canColorNode(subgraph, cNodes, edge->n1, 0) &&
        canColorNode(subgraph, cNodes, edge->n2, 1)
        ) {
        // edge is being used
        
        // indication to clean after search
        int n1C = cNodes[edge->n1] == 0 ? 0 : 1;
        int n2C = cNodes[edge->n2] == 1 ? 0 : 1;

        // color the nodes
        cNodes[edge->n1] = 0;
        cNodes[edge->n2] = 1;

        // add the edge to the graph
        subgraph->addEdge(edge);

        // continue the search with this setting
        searchBiCoSubgraphs(graph, subgraph, results, cNodes, trashWeights, curEdgeId+1);

        // remove edge from graph
        subgraph->removeEdge(edge);

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
        
        // indication to clean after search
        int n1C = cNodes[edge->n1] == 1 ? 0 : 1;
        int n2C = cNodes[edge->n2] == 0 ? 0 : 1;

        // color the nodes
        cNodes[edge->n1] = 1;
        cNodes[edge->n2] = 0;

        // add the edge to the graph
        subgraph->addEdge(edge);

        // continue the search with this setting
        searchBiCoSubgraphs(graph, subgraph, results, cNodes, trashWeights, curEdgeId+1);

        // remove edge from graph
        subgraph->removeEdge(edge);

        cNodes[edge->n1] = n1C == 1 ? -1 : 1;
        cNodes[edge->n2] = n2C == 1 ? -1 : 0;
    }

    // dont use this edge
    trashWeights += edge->weight;
    searchBiCoSubgraphs(graph, subgraph, results, cNodes, trashWeights, curEdgeId+1);
    //std::cout << "out" << "\n";
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

    // get the results
    searchBiCoSubgraphs(graph, subgraph, results, cNodes, 0, 0);

    // delete unnecesary
    delete[] cNodes;
    //for (int i = 0; i < graph->getEdgesNum(); i++) delete graph->getEdges()[i];
    //delete[] graph->getEdges();
    delete subgraph;

    // return the result
    return results;
}
