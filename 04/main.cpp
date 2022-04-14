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
int main(int argc, char **argv) {
    // test if input is correct
    if (argc != 2) {
        cout << "Incorrect input, just name of graph is needed" << "\n";
        return 1;
    }

    getMaxBiparSubgraph(argc, argv);

    printf("konec\n");

    return 0;
}