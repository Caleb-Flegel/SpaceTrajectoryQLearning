#include<vector>
#include<string>
#include <random>   // for std::mt19937_64 object
#include"../Config_Constants/constants.h"
#include"../Runge_Kutta/rkParameters.h"
#include"../Motion_Eqns/elements.h"

State::State() {
 //Set all increment values in the stateParams array to 0
    for (int i = 0; i < OPTIM_VARS; i++){
        stateParams[i] = 0;
    }
}

State::State(const cudaConstants * cConstants, std::mt19937_64 & rng) {
    //Randomly set all increment values in the stateParams
    for (int i = 0; i < OPTIM_VARS; i++){
        stateParams[i] = rng() % cConstants->num_increments;
    }
}

// Returns a vector of strings saying which actions the spacecraft could take
std::vector<std::string> State::getPossibleActions(const cudaConstants * cConstants){
    //a vector to hold the legal changes one can make to the increment value for the spacecraft
    std::vector<std::string> legalActions;
    
    //Go thru all indexes of the stateParams array and check to see if they are within the range of an increment/decrement action
    //  increment val >= 1 to allow for decrement
    //  increment val <  cConstants->num_increments to allow for increasing
    for (int i = 0; i < OPTIM_VARS; i++){
        //Checks if gamma 0 can be decremented or incremented
        if (i == GAMMA0_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA0");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA0");
            }
        }
        //Checks if gamma 1 can be decremented or incremented
        else if (i == GAMMA1_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA1");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA1");
            }
        }
        //Checks if gamma 2 can be decremented or incremented
        else if (i == GAMMA2_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA2");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA2");
            }
        }
        //Checks if gamma 3 can be decremented or incremented
        else if (i == GAMMA3_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA3");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA3");
            }
        }
        //Checks if gamma 4 can be decremented or incremented
        else if (i == GAMMA4_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA4");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA4");
            }
        }
        //Checks if gamma 5 can be decremented or incremented
        else if (i == GAMMA5_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA5");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA5");
            }
        }
        //Checks if gamma 6 can be decremented or incremented
        else if (i == GAMMA6_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_GAMMA6");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_GAMMA6");
            }
        }
        //Checks if tau 0 can be decremented or incremented
        else if (i == TAU0_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_TAU0");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_TAU0");
            }
        }
        //Checks if tau 1 can be decremented or incremented
        else if (i == TAU1_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_TAU1");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_TAU1");
            }
        }
        //Checks if tau 2 can be decremented or incremented
        else if (i == TAU2_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_TAU2");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_TAU2");
            }
        }
        //Checks if alpha can be decremented or incremented
        else if (i == ALPHA_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_ALPHA");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_ALPHA");
            }
        }
        //Checks if beta can be decremented or incremented
        else if (i == BETA_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_BETA");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_BETA");
            }
        }
        //Checks if zeta can be decremented or incremented
        else if (i == ZETA_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_ZETA");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_ZETA");
            }
        }  
        //Checks if trip time can be decremented or incremented 
        else if (i == TRIPTIME_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_TRIPTIME");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_TRIPTIME");
            }
        }
        //Checks if coast 0 can be decremented or incremented
        else if (i == COAST0_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_COAST0");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_COAST0");
            }
        }
        //Checks if coast 1 can be decremented or incremented
        else if (i == COAST1_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_COAST1");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_COAST1");
            }
        }
        //Checks if coast 2 can be decremented or incremented
        else if (i == COAST2_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_COAST2");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_COAST2");
            }
        }
        //Checks if coast 3 can be decremented or incremented
        else if (i == COAST3_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_COAST3");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_COAST3");
            }
        }
        //Checks if coast 4 can be decremented or incremented
        else if (i == COAST4_OFFSET){
            if (stateParams[i] >= 1){
                legalActions.push_back("decrement_COAST4");
            }
            if (stateParams[i] < cConstants->num_increments){
                legalActions.push_back("increment_COAST4");
            }
        }
        /*else {
            std::cout << "Something weird is going on..." << std::endl;
        }*/
    }

    //return the vector of strings containing the possible actions the spacecraft can chose from 
    return legalActions;
}

//Calculates a new state based on a (non simulate) selected action
std::array<int, OPTIM_VARS> State::getNewState(const std::string& action){
    //Based on the action increment or decrement the respective var increment value
    //  update class's stateParams array
    std::array<int, OPTIM_VARS> newParams;
    newParams = stateParams;

    if (action == "decrement_GAMMA0"){
        newParams[GAMMA0_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA0"){
        newParams[GAMMA0_OFFSET] += 1;
    }
    else if (action == "decrement_GAMMA1"){
        newParams[GAMMA1_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA1"){
        newParams[GAMMA1_OFFSET] += 1;
    }
    else if (action == "decrement_GAMMA2"){
        newParams[GAMMA2_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA2"){
        newParams[GAMMA2_OFFSET] += 1;
    }
    else if (action == "decrement_GAMMA3"){
        newParams[GAMMA3_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA3"){
        newParams[GAMMA3_OFFSET] += 1;
    }
    else if (action == "decrement_GAMMA4"){
        newParams[GAMMA4_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA4"){
        newParams[GAMMA4_OFFSET] += 1;
    }
    else if (action == "decrement_GAMMA5"){
        newParams[GAMMA5_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA5"){
        newParams[GAMMA5_OFFSET] += 1;
    }
    else if (action == "decrement_GAMMA6"){
        newParams[GAMMA6_OFFSET] -= 1;
    }
    else if (action == "increment_GAMMA6"){
        newParams[GAMMA6_OFFSET] += 1;
    }
    else if (action == "decrement_TAU0"){
        newParams[TAU0_OFFSET] -= 1;
    }
    else if (action == "increment_TAU0"){
        newParams[TAU0_OFFSET] += 1;
    }
    else if (action == "decrement_TAU1"){
        newParams[TAU1_OFFSET] -= 1;
    }
    else if (action == "increment_TAU1"){
        newParams[TAU1_OFFSET] += 1;
    }
    else if (action == "decrement_TAU2"){
        newParams[TAU2_OFFSET] -= 1;
    }
    else if (action == "increment_TAU2"){
        newParams[TAU2_OFFSET] += 1;
    }
    else if (action == "decrement_ALPHA"){
        newParams[ALPHA_OFFSET] -= 1;
    }
    else if (action == "increment_ALPHA"){
        newParams[ALPHA_OFFSET] += 1;
    }
    else if (action == "decrement_BETA"){
        newParams[BETA_OFFSET] -= 1;
    }
    else if (action == "increment_BETA"){
        newParams[BETA_OFFSET] += 1;
    }
    else if (action == "decrement_ZETA"){
        newParams[ZETA_OFFSET] -= 1;
    }
    else if (action == "increment_ZETA"){
        newParams[ZETA_OFFSET] += 1;
    }
    else if (action == "decrement_TRIPTIME"){
        newParams[TRIPTIME_OFFSET] -= 1;
    }
    else if (action == "increment_TRIPTIME"){
        newParams[TRIPTIME_OFFSET] += 1;
    }
    else if (action == "decrement_ZETA"){
        newParams[ZETA_OFFSET] -= 1;
    }
    else if (action == "increment_ZETA"){
        newParams[ZETA_OFFSET] += 1;
    }
    else if (action == "decrement_COAST0"){
        newParams[COAST0_OFFSET] -= 1;
    }
    else if (action == "increment_COAST0"){
        newParams[COAST0_OFFSET] += 1;
    }
    else if (action == "decrement_COAST1"){
        newParams[COAST1_OFFSET] -= 1;
    }
    else if (action == "increment_COAST1"){
        newParams[COAST1_OFFSET] += 1;
    }
    else if (action == "decrement_COAST2"){
        newParams[COAST2_OFFSET] -= 1;
    }
    else if (action == "increment_COAST2"){
        newParams[COAST2_OFFSET] += 1;
    }
    else if (action == "decrement_COAST3"){
        newParams[COAST3_OFFSET] -= 1;
    }
    else if (action == "increment_COAST3"){
        newParams[COAST3_OFFSET] += 1;
    }
    else if (action == "decrement_COAST4"){
        newParams[COAST4_OFFSET] -= 1;
    }
    else if (action == "increment_COAST4"){
        newParams[COAST4_OFFSET] += 1;
    }
    
    return newParams;
}

//Calculates simulation params based on the current state increments
rkParameters<double> State::getSimVal(const cudaConstants * cConstants){
    //Create new rkParams var
    //Use cConstants->minSimVals, cConstants->maxSimVals, and cConstants->numIncrements to calculate the value of each increment 
    //Multiply the increment value by the number of increments for the variable to calculate the actual value of the state 
    //Do this for 19 optim vars
    //use these to calculate rkParams.y0 (see child.cpp constructor)
    rkParameters<double> calcParams();

    //Var to store the increment value     
    double incVal = 0;

    //Go through each optimVar to calculate it's value
    //Get gamma incVal
    incVal = (cConstants->maxSimVals[GAMMA_OFFSET] - cConstants->maxSimVals[GAMMA_OFFSET])/cConstants->num_increments;

    //Calc gamma 0
    calcParams.coeff.gamma[0] = stateParams[GAMMA0_OFFSET]*incVal;
    //Calc gamma 1
    calcParams.coeff.gamma[1] = stateParams[GAMMA1_OFFSET]*incVal;
    //Calc gamma 2
    calcParams.coeff.gamma[2] = stateParams[GAMMA2_OFFSET]*incVal;
    //Calc gamma 3
    calcParams.coeff.gamma[3] = stateParams[GAMMA3_OFFSET]*incVal;
    //Calc gamma 4
    calcParams.coeff.gamma[4] = stateParams[GAMMA4_OFFSET]*incVal;
    //Calc gamma 5
    calcParams.coeff.gamma[5] = stateParams[GAMMA5_OFFSET]*incVal;
    //Calc gamma 6
    calcParams.coeff.gamma[6] = stateParams[GAMMA6_OFFSET]*incVal;


    //Get tau incVal
    incVal = (cConstants->maxSimVals[TAU_OFFSET] - cConstants->maxSimVals[TAU_OFFSET])/cConstants->num_increments;

    //Calc tau 0
    calcParams.coeff.tau[0] = stateParams[TAU0_OFFSET]*incVal;
    //Calc tau 1
    calcParams.coeff.tau[1] = stateParams[TAU1_OFFSET]*incVal;
    //Calc tau 2
    calcParams.coeff.tau[2] = stateParams[TAU2_OFFSET]*incVal;


    //Get alpha incVal
    incVal = (cConstants->maxSimVals[ALPHA_OFFSET] - cConstants->maxSimVals[ALPHA_OFFSET])/cConstants->num_increments;
    //Calc alpha val
    calcParams.alpha = stateParams[ALPHA_OFFSET]*incVal;
    
    //Get beta incVal
    incVal = (cConstants->maxSimVals[BETA_OFFSET] - cConstants->maxSimVals[BETA_OFFSET])/cConstants->num_increments;
    //Calc beta val
    calcParams.beta = stateParams[BETA_OFFSET]*incVal;
    
    //Get zeta incVal
    incVal = (cConstants->maxSimVals[ZETA_OFFSET] - cConstants->maxSimVals[ZETA_OFFSET])/cConstants->num_increments;
    //Calc zeta val
    calcParams.zeta = stateParams[ZETA_OFFSET]*incVal;

    //Get triptime incVal
    incVal = (cConstants->maxSimVals[TRIPTIME_OFFSET] - cConstants->maxSimVals[TRIPTIME_OFFSET])/cConstants->num_increments;
    //Calc triptime val
    calcParams.tripTime = stateParams[TRIPTIME_OFFSET]*incVal;

    
    //Get coast incVal
    incVal = (cConstants->maxSimVals[COAST_OFFSET] - cConstants->maxSimVals[COAST_OFFSET])/cConstants->num_increments;

    //Calc gamma 0
    calcParams.coeff.coast[0] = stateParams[COAST0_OFFSET]*incVal;
    //Calc gamma 1
    calcParams.coeff.coast[1] = stateParams[COAST1_OFFSET]*incVal;
    //Calc gamma 2
    calcParams.coeff.coast[2] = stateParams[COAST2_OFFSET]*incVal;
    //Calc gamma 3
    calcParams.coeff.coast[3] = stateParams[COAST3_OFFSET]*incVal;
    //Calc gamma 4
    calcParams.coeff.coast[4] = stateParams[COAST4_OFFSET]*incVal;

    elements<double> earth = launchCon->getCondition(calcParams.tripTime); //get Earth's position and velocity at launch

    calcParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
        earth.r+ESOI*cos(calcParams.alpha),
        earth.theta+asin(sin(M_PI-calcParams.alpha)*ESOI/earth.r),
        earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
        earth.vr+cos(calcParams.zeta)*sin(calcParams.beta)*cConstants->v_escape, 
        earth.vtheta+cos(calcParams.zeta)*cos(calcParams.beta)*cConstants->v_escape,
        earth.vz+sin(calcParams.zeta)*cConstants->v_escape);

    return calcParams;
};