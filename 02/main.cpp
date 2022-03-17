/**
 * @file main.cpp
 * @author David Omrai (omraidav@fit.cvut.cz)
 * @brief First ni-pdp project sequence algorithm
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <string>
#include "bipartSubSearcher.h"

using namespace std;

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


    Graph* graph = new Graph(graphName);

    if (graph->getMatrix() == NULL) {
        delete graph;
        return 1;
    }
   

    Results* results = getMaxBiparSubgraph(graph);

    if (results->results != NULL) {
        cout << "RESULTS" << endl;
        
        ResultNode* tmpRes = results->results;
        while (tmpRes != NULL) {
            cout << "Weight sum: " << endl;
            cout << tmpRes->weight << endl;
            cout << "Node colors:" << endl;
            for (int i = 0; i < graph->getNodesNum(); i++) {
                cout << i << ":" << tmpRes->cNodes[i] << ", ";
            }
            cout << '\n';

            // cout << "Edges:" << endl;
            // for (int i = 0; i < tmpRes->graph->getNodesNum(); i++) {
            //     for (int j = 0; j < tmpRes->graph->getNodesNum(); j++) {
            //         cout << " | " << setw(4) << tmpRes->graph->getEdgeWeight(i, j);
            //     }
            //     cout << " |" << endl;
            // }

            cout << "Edges:" << "\n";
            for (int i = 0; i < graph->getNodesNum(); i++) {
                for (int j = 0; j < graph->getNodesNum(); j++) {
                    cout << "  |   " << setw(4) << tmpRes->graphWeights[i*graph->getNodesNum()+j];
                }
                cout << "  |" << "\n";
            }

            tmpRes = tmpRes->next;
            delete results->results;
            results->results = tmpRes;

            cout << "-----------------------------------------------------------" << '\n';
        }
    }

    delete results;
    delete graph;

    return 0;
}