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
        #pragma omp critical 
        {
            // Is new solution better
            if ((results->results == NULL && isConnected(subgraph) == true) || 
                (results->results != NULL && results->results->weight <= subgraph->getWeightsSum() && isConnected(subgraph) == true)) {
                addResult(results, subgraph, cNodes, graph->getEdgesNum());
            };
        }
        delete subgraph;
        delete[] cNodes;
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
            #pragma omp task
            {
                searchColorCombination(graph, newSubgraph, results, newCNodes, trashWeights, curEdgeId+1, edge, 0, 1);
            }
            //#pragma omp taskwait
    }

    // if can be colored edge 1 0 and delete edge after you are done
    if (
        cNodes[edge->n1] != 0 && cNodes[edge->n2] != 1 &&
        canColorNode(subgraph, cNodes, edge->n1, 1) &&
        canColorNode(subgraph, cNodes, edge->n2, 0)
        ) {
            Graph* newSubgraph = subgraph->copyGraph();
            int* newCNodes = copyVector(cNodes, graph->getNodesNum());
            #pragma omp task
            {
                searchColorCombination(graph, newSubgraph, results, newCNodes, trashWeights, curEdgeId+1, edge, 1, 0);
            }
            //#pragma omp taskwait
    }

    // dont use this edge
    trashWeights += edge->weight;
    Graph* newSubgraph = subgraph->copyGraph();
    int* newCNodes = copyVector(cNodes, graph->getNodesNum());
    #pragma omp task
    {
        searchBiCoSubgraphs(graph, newSubgraph, results, newCNodes, trashWeights, curEdgeId+1);
    }
    //#pragma omp taskwait
    // // clean the mess
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
    // omp_set_num_threads(4);
    #pragma omp parallel shared(graph, results)
    {
        // get the results
        #pragma omp single
        {
            int ID = omp_get_thread_num();
            std::cout << "NUM THREADS: " << ID << "\n";
            searchBiCoSubgraphs(graph, subgraph, results, cNodes, 0, 0);
        }
    }
    
    // return the result
    return results;
}
