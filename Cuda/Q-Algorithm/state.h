#ifndef STATE_H
#define STATE_H

#include<vector>
#include<string>
#include<map>

#include "../Config_Constants/constants.h"
//#include "../Q-Algorithm/state.cpp"

// Class to handle the state
//  Includes storing the increment array, calculating the start params,
//  and finding next state actions
struct State {
    //Stores the base state increments
    // int initStateParams[OPTIM_VARS];

    // Array to store the increments for each start param
    int stateParams[OPTIM_VARS];

    //Map to store a state's value
    std::map<int[OPTIM_VARS], double> values;

    //Base constructor
    State();

    //Construction 
    State(const cudaConstants * cConstants, std::mt19937_64 & rng);

    // Returns a vector of strings saying which actions the spacecraft could take
    std::vector<std::string> getPossibleActions(const cudaConstants * cConstants);

    //Calculates a new state based on a (non simulate) selected action
    std::array<int, OPTIM_VARS> getNewState(const std::string& action);

    //Calculates simulation params based on the current state increments
    rkParameters<double> getSimVal(const cudaConstants * cConstants);
};

#include "state.cpp"

#endif