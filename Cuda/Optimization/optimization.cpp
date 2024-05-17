// Spaceflight Optimization Project using CUDA and a genetic algorithm

#include "../Planet_calculations/planetInfo.h"  // For launchCon and EarthInfo(); includes elements, runge_kutta, & config.h
//#include "../Q-Algorithm/adult.h" // For adult structs, paths to rkParameters for randomParameters(); includes child, rkParameters, & planetInfo
#include "../Q-Algorithm/child.h" // For child structs, paths to rkParameters for randomParameters(); includes config.h, rkParameters, & planetInfo
#include "../Output_Funcs/output.h" // For terminalDisplay(), recordGenerationPerformance(), and finalRecord()
#include "../Runge_Kutta/runge_kutta.h" // for testing rk4simple; includes calcFourier, motion_eqns, child, & gpuMem
//#include "../Q-Algorithm/ga_crossover.h" // for selectSurvivors() and newGeneration(); includes constants.h, adult, & child
//#include "../Q-Algorithm/genetic_algorithm.h" //For functions that set up new generations; includes constants.h, adult, child, ga_crossover, & sort
//#include "../Q-Algorithm/referencePoints.h" //For the ReferencePoints class which deals with calculating reference points and setting adult's rarity
//#include "../Q-Algorithm/sort.h" //For functions that will allow for sorting of the adult arrays by giving them ranks and distances; includes constants.h & adult
//#include "../Q-Algorithm/anneal.h" //For all the annealing functions; includes constants.h & adult
//#include "../Runge_Kutta/gpuMem.cuh" // for initializing and deallocating; includes child, rkParameters, & config.h

#include <iostream> // cout
#include <iomanip>  // used for setw(), sets spaces between values output
#include <random>   // for std::mt19937_64 object
#include <vector>   // allows us to use vectors instead of just arrays
#include <array>
#include <string>
#include <chrono>
#include <algorithm>


//----------------------------------------------------------------------------------------------------------------------------
// ** Assumes pool is sorted array of Adults **
// Used in determining if main optimize loop continues
// Input: oldAdults - this generation of Adults, defined/initialized in optimimize
//        cConstants - struct holding config values, used for accessing best_count value and objectives list
// Output: Returns true if top best_count adults within the pool are within the tolerance
bool checkTolerance(const Child& ind, const cudaConstants* cConstants);

//----------------------------------------------------------------------------------------------------------------------------
// Adjusts epsilon based on how many generations have passed and sets the 
// void resetAndAdjustEps(const int generation, State& bestState, State& curState);

//Calculates the new epsilon value
double newEpsilon(const cudaConstants* cConstants, const int& generation);

//----------------------------------------------------------------------------------------------------------------------------
// TEST / LIKELY TEMPORARY FUNCTION
// This function will find the minimum, maximum, and average distance, average pos and speed diffs, the number of duplicates,
//                         the avg age, and the avg and max birthdays of the individuals in allAdults, which will then be used for reporting
// 
// Inputs:  allAdults - array of adults that will be considered
//          objectives - the vector of this run's objectives
//          objectiveAvgValues - a vector which will hold the generations's average parameter values for each of the objectives
//          duplicateNum - the number of duplicate adults found
//          minDist - minimum distance that will be calculated
//          avgDist - the average distance that will be calculated
//          maxDist - the maximum distance that will be calculated
//          generation - the current generation
//          avgAge  - the avegrage age of the adults, relative to the current generation
//          avgBirthday - average birth generation for the adults
//          oldestBirthday - the oldest adult's birth generation
// Outputs: The arguments will be filled in with the up-to-date values for this generation
//void calculateGenerationValues (std::vector<Adult> & allAdults, const std::vector<objective> & objectives, std::vector<double> & objectiveAvgValues, int & totAssoc, int & duplicateNum, double & minDist, double & avgDist, double & maxDist, int & minSteps, int & avgSteps, int & maxSteps, const int & generation, double & avgAge, double & avgBirthday, int & oldestBirthday);

//----------------------------------------------------------------------------------------------------------------------------
// Main processing function for Genetic Algorithm
// - manages memory needs for genetic algorithm
// - deals with processing calls to CUDA callRK
// - exits when individuals converge on tolerance defined in Constants
double optimize(const cudaConstants* cConstants);

//-----------------------------------------------------------------------------------------------------------------------------
int main () {
    
    // display GPU properties and ensure we are using the right one
    // cudaDeviceProp prop;
    // cudaGetDeviceProperties(&prop, 0);
    // std::cout << "\n\nDevice Number: 0 \n";
    // std::cout << "- Device name: " << prop.name << std::endl << std::endl;
    // cudaSetDevice(0);

    // Declare the genetic constants used, with file path being used to receive initial values
    cudaConstants * cConstants = new cudaConstants(); 

    //preallocates all the memory for the varaibles used by the GPU
    //also allows the GPU to access the marsLaunchCon without reloading it everytime
    // GPUMem gpuValues;

    //Tell user when the referennce points are being created
    // std::cout << "\nCreating reference points...\n";

    //Creates the reference points for the rest of the program
    // ReferencePoints refPoints(cConstants);

    //Display the number of reference points created as a sanity check
    // std::cout << "\n" << refPoints.points.size() << " reference points created.\n";

    // Sets run0 seed, used to change seed between runs
    // Seed is set in cudaConstants: current time or passed in via config
    double zero_seed = cConstants->time_seed;
    // Perform the optimization with optimize function
    for (int run = 0; run < cConstants->run_count; run++) {
        // Adjust the time_seed so it is unique based on each run
        cConstants->time_seed = zero_seed + run*100;

        // Display contents of cConstants being used for this run and how many runs
        std::cout << *cConstants;
        std::cout << "\tPerforming run #" << run+1 << "\n\n";

        // pre-calculate a table of Earth's and Mars' position within possible mission time range
        // defined as global variable
        // accessed on the CPU when individuals are initialized
        launchCon = new PlanetInfo(cConstants, EARTH); 
        marsLaunchCon = new PlanetInfo(cConstants, MARS);
        //This ensures that we copy the correct size of marsCon to the GPU
        int marsConSize = getPlanetSize(cConstants);
        //initialize all values needed for GPU calculations
        // gpuValues.initialize(cConstants, marsConSize, marsLaunchCon->getAllPositions());

        // Call optimize with the current parameters in cConstants
        optimize(cConstants);
        
        // Deallocate launchCon info for this run as it may be using a different time range in the next run
        delete launchCon; 
        delete marsLaunchCon;
        // gpuValues.free();
    }
    // Now that the optimize function is done (assumed that optimize() also records it), deallocate memory of the cudaConstants
    delete cConstants;
    
    return 0;
}

//Returns true if top best_count adults within the oldAdults vector are within the tolerance
bool checkTolerance(const Child& ind, const cudaConstants* cConstants) {
    //Sort the vector by progress to make sure the program checks the correct adult
    // std::sort(oldAdults.begin(), oldAdults.end(), bestProgress); 

    //The function needs to check if the best adult meets the convergence tolerance for each objective
    //Iterate through the objectives
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //Check to see if the top best_count adults have met convergence for this parameter
        // for (int j = 0; j < cConstants->best_count; j++) {

            //Check to see if the adult's parameter is larger than the convergence (plus equate tolerance)
            if (ind.objTargetDiffs[i] > (cConstants->missionObjectives[i].allowedDifference + cConstants->missionObjectives[i].equateTolerance)) {
                //Return false as a parameter that needs to be minimized is larger than the convergence threshold
                return false;
            }
        // }
    }

    //If the program reaches this spot, it means all of the adult's parameters have met the convergence threshold
    //  Otherwise, the function would have already returned false
    //  Thus, the adult has converged and it is appropriate to return true
    return true; 
}

//Function that will calculate distance and birthday values for a generation
// void calculateGenerationValues (std::vector<Adult> & allAdults, const std::vector<objective> & objectives, std::vector<double> & objectiveAvgValues, int & totAssoc, int & duplicateNum, double & minDist, double & avgDist, double & maxDist, int & minSteps, int & avgSteps, int & maxSteps, const int & generation, double & avgAge, double & avgBirthday, int & oldestBirthday){
//     //Reset the dist values
//     minDist = static_cast<int>(objectives.size() * MAX_DISTANCE); //Set the min dist to the maximum possible value, so that it will be changed
//     avgDist = 0; 
//     maxDist = 0; //Set the max dist to the min possible value, so that it is garunteed to be changed

//     minSteps = 1000000;
//     avgSteps = 0;
//     maxSteps = 0;

//     //Reset the average parameter values
//     objectiveAvgValues.clear();
//     //Prep the avg value vector to have an index for each objective
//     for (int i = 0; i < objectives.size(); i++) {
//         objectiveAvgValues.push_back(0.0); 
//     }
    
//     //Reset the age values
//     avgAge = 0;
//     avgBirthday = 0; 
//     oldestBirthday = generation; //Set to generation, since it is the "newest" possible oldest birthday

//     //Find the total number reference points used
//     //Reset the total points used
//     totAssoc = 1;
//     //trackers used to determine when a new ref point is used
//     int prevPoint = 0;

//     //Sort the adults by lowest ref point index to make it easy to find the number of reference points
//     std::sort(allAdults.begin(), allAdults.end(), lowAssocPoint);
    
//     //For loop will iterate through the adult array to find the values needed
//     for (int i = 0; i < allAdults.size(); i++) {
//         //Check to see if this adult is associated with a different reference point than before
//         if (allAdults[i].associatedRefPoint != allAdults[prevPoint].associatedRefPoint) {
//             //Add one to the number of used points
//             totAssoc++; 

//             //Set the new prevPoint as this index
//             prevPoint = i;
//         }

//         //Check to see if this adult's distance is a new minimum
//         if (allAdults[i].distance < minDist) {
//             minDist = allAdults[i].distance; //Set the new min distance
//         }
//         //Check to see if this adult's distance is the new maximum
//         else if (allAdults[i].distance > maxDist) {
//             maxDist = allAdults[i].distance; //Ser the new max distance
//         }

//         //Check to see if this adult's step count is a new minimum
//         if (allAdults[i].stepCount < minSteps) {
//             minSteps = allAdults[i].stepCount; //Set the new min distance
//         }
//         //Check to see if this adult's step count is the new maximum
//         else if (allAdults[i].stepCount > maxSteps) {
//             maxSteps = allAdults[i].stepCount; //Set the new max distance
//         }

//         //Check to see if this adult's birthday is older than the current oldest
//         if (allAdults[i].birthday < oldestBirthday) {
//             oldestBirthday = allAdults[i].birthday; 
//         }
           
//         //Add to the avg distance
//         avgDist += allAdults[i].distance;
//         avgSteps += allAdults[i].stepCount;

//         //Add the adult's parameter values to the necessary spot in the objective average value vector
//         for (int j = 0; j < objectives.size(); j++) {
//             objectiveAvgValues[j] += allAdults[i].getParameters(objectives[j]);
//         }
        
//         //Add to the avg age values
//         //avgAge measures the age relative to the current generation (how many generations old), so it is generation minus the adults' birthdays 
//         avgBirthday += allAdults[i].birthday;
//         avgAge += (generation - allAdults[i].birthday);   
//     }
//     //Possible floating point roundoff error
//     //  Shouldn't matter, as we don't ever compared averages to individuals
//     //Divide the averages by the number of adults to get the averages
//     avgDist /= allAdults.size();
//     avgSteps /= allAdults.size();
//     for (int i = 0; i < objectiveAvgValues.size(); i++) {
//         objectiveAvgValues[i] /= allAdults.size();
//     }
//     avgAge /= allAdults.size(); 
//     avgBirthday /= allAdults.size(); 
    
// }

//Calculates the new epsilon value
double newEpsilon(const cudaConstants* cConstants, const int& generation) {
    // Get run progress fraction
    double runProg = generation/cConstants->max_generations;

    // Calculate the new epsilon
    double epsilonAddVal = abs(cConstants->epsilon_Final - cConstants->epsilon_Initial) * runProg;
    return cConstants->epsilon_Final - epsilonAddVal;
}

//Function that gets a new state for the individual
//  Compares all possible actions and randomly chooses either a random action or the best action and returns the new state based on that action
std::vector<int> getNextState(const cudaConstants* cConstants, Child & ind, std::mt19937_64 & rng, const double& curEpsilon) {
    //Get possible actions
    std::vector<std::string> possActions = ind.curState.getPossibleActions(cConstants);

    // std::cout << "\n" << possActions.size() << "\n";

    //Vector that will store the next states based on the actions
    std::vector<std::vector<int>> nextStates;

    //Fill the nextstates vector based on the possible actions
    for (int i = 0; i < possActions.size(); i++){
        std::vector<int> nextState = ind.curState.getNewState(possActions[i]);
        nextStates.push_back(nextState);

        // Check to see if this state has been tested yet
        if (ind.curState.values[nextState] == 0) {
            // Set to the initial value
            ind.curState.values[nextState] = cConstants->initialValue;
        } 
    }

    //Array which will hold the output state
    std::vector<int> outputState;

    //Generate random value to see if the selected next state is random or the max value
    float stateType = rng() / rng.max();

    //See if the next state should be a random or best state
    if (stateType > curEpsilon){
        //Get the max value state
        //index with the max nextState val
        int max = 0;
        //For loop to compare the rest of the next state values
        for (int i = 1; i < nextStates.size(); i++) {
            if (ind.curState.values[nextStates[i]] > ind.curState.values[nextStates[max]]) {
                max = i;
            }
        }

        //Assign the max value state
        outputState = nextStates[max];
    }
    else {
        //Get a random state
        float randState = rng() % nextStates.size();

        //Assign the random state
        outputState = nextStates[randState];
    }

    return outputState;
    
}

//update value for current state, get the next state --> loop this --> optimization loop, checking tolerance

//we have simulated the new state, the spacecraft has moved to the new state, and we update the q-value of the previous state
//prevState is basically the key for the previous state
void update_q_values(const cudaConstants* cConstants, Child & ind, const std::vector<int>& prevState){
    //we want the reward to be how close to the target the individual is
    double reward = ind.progress + cConstants->livingReward;

    //sample = reward + self.discount *self.getValue(nextState);

    //takes the known reward and any known information about the value of the function at the next state to update the value at the current state
    double sample = reward + cConstants->gamma*ind.curState.values[ind.curState.stateParams]; //is max actually defined 

    //self.Q_values_curr[(state, action)] = (1-self.alpha)*self.getQValue(state, action) + self.alpha*sample

    //then find the new Q estimate from the current one 
    ind.curState.values[prevState] = (1 - cConstants->alpha)*ind.curState.values[prevState] + cConstants->alpha*sample;
}

//Function that will facilitate the process of finding an optimal flight path
double optimize(const cudaConstants* cConstants) {
    // Not used, previously used for reporting computational performance
    double calcPerS = 0;

    time_t timeSeed = cConstants->time_seed;
    std::mt19937_64 rng(timeSeed); // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    
    std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;

    //Initialize the output object with the base folder location (..\Output_Files\)
    output genOutputs(cConstants);

    // Initial genetic anneal scalar
    //double currentAnneal = cConstants->anneal_initial;
    
    // Main set of parameters for Genetic Algorithm
    // contains all thread unique input parameters
    // The "children pool" of the current genertation
    //std::vector<Adult> newAdults; 

    //Input parameters from the previous generation. Mixes with new children to determine new inputParameters
    // The "potential parent pool" of the current generation
    // DISCLAIMER - this is mentioned as the "parent pool", but it's actually the *POTENTIAL* parent pool. The survivor pool is what's used to generate the next generation. Survivors are taken from this, so it's more accurate to call this the "Potential Parent Pool"
    //std::vector<Adult> oldAdults; 

    //the set of all old and new individuals
    //std::vector<Adult> allAdults;

    // number of current generation
    int generation = 0;

    genOutputs.recordMarsData(cConstants,generation);
    genOutputs.recordEarthData(cConstants,generation);
    // genOutputs.recordReferencePoints(cConstants, refPoints);

    // Flag for finishing the genetic process
    // set by checkTolerance()
    bool convergence = false;

    //number of errors, specifically used for diagnostic recording
    //  couting all adults in the generation - includes oldAdults and newAdults that have nan values
    // int numErrors = 0;
    // int marsErrors = 0;

    //Number of different reference points that have adults associated with them every generation
    // int totAssoc = 0;

    //Initialize variables needed for distance, number of duplicate adults, and birthday reporting
    //int duplicateNum = 0;
    //double maxDistance, minDistance, avgDistance, avgAge, avgBirthday;
    //int minSteps, avgSteps, maxSteps, oldestBirthday;

    //Inititalize variables for storing the average time per generation
    std::chrono::time_point<std::chrono::system_clock> runStartTime = std::chrono::system_clock::now();
    std::chrono::duration<float> totRunTime;

    //Vector used to report the average parameter value for each objective
    // std::vector<double> objectiveAvgValues; 

    //Creates the individuals needed for the 0th generation
    //Need to make children, then callRK, then make into adults (not currently doing that)
    //oldAdults goes into this function empty and comes out filled with num_individuals Adults
    //      these adults are either randomly generated or pulled from a file
    // createFirstGeneration(oldAdults, cConstants, rng, generation, gpuValues, refPoints); 

    //Cur epsilon value
    double curEpsilon = .5;

    //TODO (DONE): Generate initial state here
    //Generate initial state
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE INIT State-_-_-_-_-_-_-_-_-_\n\n";
    State initState(cConstants, rng);
    //Generate individual
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE INIT CHILD-_-_-_-_-_-_-_-_-_\n\n";
    Child individual(initState, cConstants);
    Child bestIndividual(initState, cConstants);



    // main gentic algorithm loop
    // - continues until checkTolerance returns true (specific number of individuals are within threshold)
    do {

        //Get new epsilon
        curEpsilon = newEpsilon(cConstants, generation);

        // std::cout << "\nINCOMING STATE:\n";
        // for (int i = 0; i < individual.curState.stateParams.size(); i++){
        //     std::cout << individual.curState.stateParams[i] << "\t";
        // }

        // Genetic Crossover and mutation occur here
        //takes in oldAdults (the potential parents) and fills newAdults with descendants of the old adults
        //oldAdults is filled with the potential parents for a generation (num_individuals size) 
        //      after the run, oldAdults should remain the same
        //newAdults is empty and will be filled with the "grown" children generated in this method
        // std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE NEW GEN-_-_-_-_-_-_-_-_-_\n\n";
        // newGeneration(oldAdults, newAdults, currentAnneal, generation, rng, cConstants, gpuValues);

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE New State-_-_-_-_-_-_-_-_-_\n\n";
        //TODO: use function to get possible actions then get the value of the adjacent actions
        std::vector<int> newState = getNextState(cConstants, individual, rng, curEpsilon);

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE Assign New State-_-_-_-_-_-_-_-_-_\n\n";
        //Save the previous state and update the child's state with the new state
        std::vector<int> prevState = individual.curState.stateParams;
        individual.curState.stateParams = newState;

        // std::cout << "\nNEW STATE:\n";
        // for (int i = 0; i < individual.curState.stateParams.size(); i++){
        //     std::cout << individual.curState.stateParams[i] << "\t";
        // }

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE Get Sim VAL-_-_-_-_-_-_-_-_-_\n\n";
        //elements<double> earth = launchCon->getCondition(calcParams.tripTime); //get Earth's position and velocity at launch
        individual.curParams = individual.curState.getSimVal(cConstants);

        elements<double> earth = launchCon->getCondition(individual.curParams.tripTime); //get Earth's position and velocity at launch

        individual.curParams.y0 = elements<double>( // calculate the starting position and velocity of the spacecraft from Earth's position and velocity and spacecraft launch angles
            earth.r+ESOI*cos(individual.curParams.alpha),
            earth.theta+asin(sin(M_PI-individual.curParams.alpha)*ESOI/earth.r),
            earth.z, // The spacecraft Individual is set to always be in-plane (no initial Z offset relative to earth) 
            earth.vr+cos(individual.curParams.zeta)*sin(individual.curParams.beta)*cConstants->v_escape, 
            earth.vtheta+cos(individual.curParams.zeta)*cos(individual.curParams.beta)*cConstants->v_escape,
            earth.vz+sin(individual.curParams.zeta)*cConstants->v_escape);

        // std::cout<< "\n" << individual.curParams << "\n";

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE callRK-_-_-_-_-_-_-_-_-_\n\n";
        //TODO: sim here to get progress
        double timeInitial = 0;
        callRKBasic(individual, cConstants->rk_tol, cConstants, marsLaunchCon->getAllPositions(), timeInitial);
        
        
        // td::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE update_q_values-_-_-_-_-_-_-_-_-_\n\n";

        //TODO: (possibly in the last func) choose either the best or a random next state
        update_q_values(cConstants, individual, prevState);
        
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE PREP PARENTS-_-_-_-_-_-_-_-_-_\n\n";
        //fill oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected (numErrors is for all adults in the generation - the oldAdults and the newAdults)
        //allAdults will be filled with the last generation's set of parents and their offspring when this starts (sorted by rankDistanceSort)
        //      by the end of this function, it will be filled with the new generation's set of parents and children sorted by rankDistanceSort
        //newAdults goes in with the "grown" children created in new generation (size num_individuals)
        //      by the end of the function, it is cleared
        //oldAdults goes in with the pool of potential parents that may have generated the newAdults
        //      by the end of the function, it is filled with the best num_individuals adults from allAdults (sorted by rankDistanceSort) 
        //preparePotentialParents(allAdults, newAdults, oldAdults, numErrors, duplicateNum, cConstants, rng, refPoints, generation, currentAnneal, marsErrors);


        // Display a '.' to the terminal to show that a generation has been performed
        // This also serves to visually seperate the terminalDisplay() calls across generations 
        std::cout << '.';

        //TODO: Add function to adjust alpha, gamma, epsilon levels 

        //Perform utitlity tasks (adjusting anneal and reporting data)
        //Calculate variables for birthdays and distances
        //calculateGenerationValues(allAdults, cConstants->missionObjectives, objectiveAvgValues, totAssoc, duplicateNum, minDistance, avgDistance, maxDistance, minSteps, avgSteps, maxSteps, generation, avgAge, avgBirthday, oldestBirthday);

        //reset sort from calculateGenerationValues function
        //mainSort(allAdults, cConstants, allAdults.size());

        //Assumes oldAdults is in rankDistance order
        //changeAnneal (oldAdults, cConstants, currentAnneal, generation);

        //Get the run's total time
        totRunTime = (std::chrono::system_clock::now() - runStartTime);

        if (individual.progress >= bestIndividual.progress){
            bestIndividual = individual;
        }

        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE RECORD-_-_-_-_-_-_-_-_-_\n\n";
        //Print out necessary info for this generation
        genOutputs.printGeneration(cConstants, individual, bestIndividual, generation, (totRunTime.count()/(generation+1)));

        // Before replacing new adults, determine whether all are within tolerance
        // Determines when loop is finished
        //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE CONVERGENCE CHECK-_-_-_-_-_-_-_-_-_\n\n";
        convergence = checkTolerance(individual, cConstants); //QL: done

        //Reset the progress sort from checkTolerance
        //mainSort(oldAdults, cConstants, oldAdults.size());   

        //Reset if the answer hasn't been found
        if (!convergence && (generation % cConstants->resetGenNum) == 0){
            //Set the individual's state params to the best individual's params
            std::cout << "\nState Reset!\n";
            individual.curState.stateParams = bestIndividual.curState.stateParams;
        }   
        
        //Increment the generation counter
        ++generation;

        //Reset child
        individual= Child(individual.curState, cConstants);
    
        //Loop exits based on result of checkTolerance and if max_generations has been hit
    } while ( !convergence && generation < cConstants->max_generations);

    //Handle the final printing
    genOutputs.printFinalGen(cConstants, individual, bestIndividual, convergence, generation, (totRunTime.count()/(generation+1))); 

    return calcPerS;
}

