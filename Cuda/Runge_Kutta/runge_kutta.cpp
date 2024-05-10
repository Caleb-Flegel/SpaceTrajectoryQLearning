
#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI

#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions
#include <fstream> // for outputing to .csv file
#include <chrono>
#include <algorithm> // sort(), shuffle()
#include <random>

#include "../Motion_Eqns/motion_equations.h"
#include "..\Runge_Kutta\rkParameters.h"
#include "..\Config_Constants\config.h"

/* //The following is outdated as of summer 2023, see the RK4Sim in ruge_kuttaCUDA for the used simulation function
template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y_new, 
                               const T & absTol, coefficients<T> coeff, T *gamma,  T *tau, int & lastStep, T *accel_output, T *fuelSpent, const T & wetMass, const cudaConstants* cConstant, elements<T> *marsIndex) {

    //Create a sim status which indidicates the run has not started
    SIM_STATUS simStatus = INITIAL_SIM;

    T curTime; // setting time equal to the start time
    T startTime;
    T endTime;
    T t_SOI;//Time stamp at a SOI boundary
    int n = 0; // setting the initial step number equal to 0

    T massFuelSpent;  //mass of fuel expended (kg), set to 0 initially

    elements<T> u;//Current orbital elements (r, theta, z, vr, vtheta, vz)
    elements<T> y_SOI;//Orbital elements at a SOI boundary

    //While loop will continute the calculations until the run has completed
    while (simStatus != COMPLETED_SIM) {

        std::cout << "\nNEW sim cycle with stats: " << simStatus << "\n";
        
        if (simStatus == INITIAL_SIM){//Starting from launch conditions
            //Child has not been simulated, set the initial curTime to the start time of the simulation
            curTime = timeInitial;
            //Set the start time to the total trip time
            startTime =  timeInitial;
            // start with the initial conditions of the spacecraft
            u = y0;

            //No fuel has been spent initially
            massFuelSpent = 0;
        }
        else{//Has been partially simulated (has entered a SOI), reset initial conditions
            curTime = t_SOI;
            //Set the start time to this simulation's start time
            startTime = t_SOI;
            // start with the initial conditions of the spacecraft
            u = y_SOI;
        }
        //Set the proper endTime
        //If it is an SOI run, the endTime will be curTime plus the estimated assist/orbital time
        //If not, end time will equal timeFinal
        if (simStatus == INSIDE_SOI) {
            //If this is a mission inside a SOI, set the final to be a gravAssistTime difference from the start time
            endTime = startTime + cConstant->gravAssistTime;
            //endTime = curTime + ((timeFinal - startTime) * cConstant->gravAssistTimeFrac);

            //Check to make sure that endTime is not further than tripTime
            if(endTime > timeFinal) {
                endTime = timeFinal;
            }
        }
        //In all other scenerios, the end time is the triptime
        else {
            endTime = timeFinal;
        }

        thruster<T> thrust(cConstant);

        bool coast; //=1 means thrusting, from calc_coast()

        elements<T> error; // error needs to be defined, used in calc_scalingFactor

        while (curTime < endTime) { // iterate until time is equal to the stop time
            times[n] = curTime;

            // Check the thruster type before performing calculations, at time 0
            if (cConstant->thruster_type == thruster<double>::NO_THRUST) {
                // gamma[n] = tau[n] = accel_output[n] = fuelSpent[n] = 0;
                coast = accel_output[n] = 0;
            }
            else {
                // determines whether the craft is coasting
                coast = calc_coast(coeff, curTime, timeFinal, thrust);
                // array of acceleration for binary output
                accel_output[n] = calc_accel(u.r,u.z, thrust, massFuelSpent, stepSize, coast, wetMass, cConstant);
                // array of gamma for binary output
                gamma[n] = calc_gamma(coeff, curTime, timeFinal);
                // array of tau for binary output
                tau[n] = calc_tau(coeff,curTime, timeFinal);
                // array of fuel spent for binary output
                fuelSpent[n] = massFuelSpent;
            }

            //WARNING/NOTE: The indexing of n (not getCondition) for the error here should be double checked
            elements<double> mars = (*marsLaunchCon).getCondition(timeFinal - curTime); //gets Mars' position relative to the Sun
            marsIndex[n] = mars;

            //calculate the distance between the spacecraft and Mars
            double marsCraftDist = sqrt(pow(mars.r, 2) + pow(u.r, 2) + pow(u.z - mars.z, 2) - (2*u.r*mars.r*cos(mars.theta-u.theta)));

            //calculate new position
            //rkCalc(curTime, timeFinal, stepSize, u, coeff, accel_output[n], error, mars, marsCraftDist);
            rkCalc(curTime, timeFinal, stepSize, u, coeff, accel_output[n], error, mars, marsCraftDist);

            //array of time output as t         
            curTime += stepSize;

            ////This is the way that stepSize was calculated in rk4CUDASim
            //stepSize *= calc_scalingFactor(u-error,error,absTol, cConstant->doublePrecThresh); // Alter the step size for the next iteration

            //// The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            //if (stepSize > (endTime - startTime) / cConstant->min_numsteps) {
            //    stepSize = (endTime - startTime) / cConstant->min_numsteps;
            //}
            //else if (stepSize < (endTime - startTime) / cConstant->max_numsteps) {
            //    stepSize = (endTime - startTime) / cConstant->max_numsteps;
            //}
            
            stepSize = (endTime - startTime) / cConstant->max_numsteps;
            
            // shorten the last step to end exactly at the end time
            if ( (curTime + stepSize) > endTime) {
                stepSize = (endTime - curTime); // shorten the last step to end exactly at time final
            }

            //Check to see if the curTime is less than endTime
            if (curTime < endTime) {
                //If so, check to see if the child triggers the conditions for a new run

                //See if child has entered MSOI
                if (marsCraftDist < MSOI*cConstant->MSOI_scale && simStatus != INSIDE_SOI) {
                    //Reset the conditions
                    t_SOI = curTime; 
                    y_SOI = u; 

                    //Set the child status to inside SOI
                    simStatus = INSIDE_SOI;                

                    std::cout << "\n Enters SOI at step" << n <<" and time" << curTime << "\n";

                    break;
                }

                //Check if child has exited MSOI after being inside it
                if (marsCraftDist > MSOI*cConstant->MSOI_scale && simStatus == INSIDE_SOI) {
                    //Reset the conditions
                    t_SOI = curTime; 
                    y_SOI = u;

                    //Set the child status to outide SOI
                    simStatus = OUTSIDE_SOI;

                    std::cout << "\n Exits SOI at step" << n <<" and time" << curTime << "\n";
                    
                    break;
                }
            } // end if (curTime < endTime)

            //The time is at or past the endTime 
            else {
                //if the endTime is at the final time, the simulation is done
                if (endTime >= timeFinal) {
                    std::cout << "\nCompleted Sim status assigned at step " << n << " at time " << curTime << "\n";
                    simStatus = COMPLETED_SIM;
                }
                //if not, it means that the simulation is in a SOI and hasn't escaped the SOI by the time the estimated orbit/assist time
                //set up the next simulation pass by resetting time initial to the current time
                else {
                    t_SOI = curTime;
                    y_SOI = u;
                }
            }//end of if-else (curTime < endTime) after updating curTime

            y_new[n] = u;              
            n++;
        } //end of while (curTime < endTime) before updating curTime
    }//end of !COMPLETED_SIM

    std::cout << "\nCOMPLETED SIM\n";
     
    lastStep = n-1;//Array index begins at 0; n was incremented after y_new was assigned.

    std::cout << "\n\tRK4SYS: " << n << "\n\n";

    // Test outputs to observe difference between rk4sys results with CUDA runge-kutta results
    if (cConstant->orbitalSpeed == NOT_APPLICABLE){//change the calculation if it is NOT (?) an orbital mission 
        std::cout << "rk4sys posDiff: " << sqrt(pow(cConstant->r_fin_target - y_new[lastStep].r, 2) + pow(cConstant->r_fin_target * cConstant->theta_fin_target - y_new[lastStep].r * fmod(y_new[lastStep].theta, 2 * M_PI), 2) + pow(cConstant->z_fin_target - y_new[lastStep].z, 2)) << std::endl;
        std::cout << "rk4sys speedDiff: " << sqrt(pow(cConstant->vr_fin_target - y_new[lastStep].vr, 2) + pow(cConstant->vtheta_fin_target - y_new[lastStep].vtheta, 2) + pow(cConstant->vz_fin_target - y_new[lastStep].vz, 2));
    }else{
            std::cout << "rk4sys orbitalposDiff: " << abs(sqrt(pow(cConstant->r_fin_target, 2) + pow(y_new[lastStep].r, 2) + pow(y_new[lastStep].z - cConstant->z_fin_target, 2) - (2*y_new[lastStep].r*cConstant->r_fin_target*cos(cConstant->theta_fin_target-y_new[lastStep].theta))) - cConstant->orbitalRadius) << std::endl;
            std::cout << "rk4sys orbitalspeedDiff: " << sqrt(abs((pow(cConstant->vr_fin_target - y_new[lastStep].vr, 2) + pow(cConstant->vtheta_fin_target - y_new[lastStep].vtheta, 2) + pow(cConstant->vz_fin_target - y_new[lastStep].vz, 2)) - pow(cConstant->orbitalSpeed, 2)));

            // std::cout << "\n\nRK4SYS r initial of craft: " << y0.r;
            // std::cout << "\nRK4SYS theta initial of craft: " << y0.theta;
            // std::cout << "\nRK4SYS z initial of craft: " << y0.z;
            //std::cout << "\n\nRK4SYS final R of Target: " << cConstant->r_fin_target;
            //std::cout << "\nRK4SYS final Theta of Target: " << cConstant->theta_fin_target;
            //std::cout << "\nRK4SYS final Z of Target: " << cConstant->z_fin_target;

            std::cout << "\n\nRK4SYS r final of craft: " << y_new[lastStep].r;
            std::cout << "\nRK4SYS theta final of craft: " << y_new[lastStep].theta;
            std::cout << "\nRK4SYS z final of craft: " << y_new[lastStep].z;

            // std::cout << "\n\nRK4SYS vr initial of craft: " << y0.vr;
            // std::cout << "\nRK4SYS vtheta initial of craft: " << y0.vtheta;
            // std::cout << "\nRK4SYS vz initial of craft: " << y0.vz;

            //std::cout << "\n\nRK4SYS final vr of Target: " << cConstant->vr_fin_target;
            //std::cout << "\nRK4SYS final vTheta of Target: " << cConstant->vtheta_fin_target;
            //std::cout << "\nRK4SYS final vz of Target: " << cConstant->vz_fin_target;
            // std::cout << "\n\nRK4SYS vr final of craft: " << y_new[lastStep].vr;
            // std::cout << "\nRK4SYS vtheta final of craft: " << y_new[lastStep].vtheta;
            // std::cout << "\nRK4SYS vz final of craft: " << y_new[lastStep].vz<<"\n";

            // std::cout << "\nRK4SYS fuel spent: " << fuelSpent[lastStep]<<"\n";

            // std::cout << "\n\nRK4SYS initial r of craft: " << y0.r;
            // std::cout << "\nRK4SYS initial theta of craft: " << y0.theta;
            // std::cout << "\nRK4SYS initial z of craft: " << y0.z;
            // std::cout << "\nRK4SYS initial vr of craft: " << y0.vr;
            // std::cout << "\nRK4SYS initial vtheta of craft: " << y0.vtheta;
            // std::cout << "\nRK4SYS initial vz of craft: " << y0.vz;
    }   
}*/

/*
// Called by optimize() in optimization.cu
//void callRK(const int numThreads, const int blockThreads, Child *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, const cudaConstants* cConstant, PlanetInfo *marsLaunchCon) {
void callRK(Child individual, double absTolInput, const cudaConstants* cConstant, elements<double> *marsLaunchCon, std::vector<double>& time_steps, std::vector<elements<double>>& y_steps, std::vector<double>& gamma_steps, std::vector<double>& tau_steps, std::vector<double>& accel_steps, std::vector<double>& fuel_steps, double timeInitial) {

    //Counters for the number of children who are finished being simulated and how many simulations the code has ran for this generation
    // int completedNum;
    // int simulationNum = 0;

    //is this bad?
    //individual.simStatus = INITIAL_SIM;
    
    double stepSize = (cConstant->triptime_min - timeInitial)/cConstant->max_numsteps;
    
    //Do while loop will keep simulating while there are still children to simulate
    do {

        //Reset the count of children who have finished their simulation
        // int completedNum = 0;

        // cudaEvent_t kernelStart, kernelEnd;
        // cudaEventCreate(&kernelStart);
        // cudaEventCreate(&kernelEnd);

        // // copy values of parameters passed from host onto device
        // cudaMemcpy(gpuValues.devGeneration, generation, numThreads * sizeof(Child), cudaMemcpyHostToDevice);
        // // cudaMemcpy(gpuValues.devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
        // // cudaMemcpy(gpuValues.devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);

        // GPU version of rk4Simple()
        // cudaEventRecord(kernelStart);
        //rk4CUDASim<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(gpuValues.devGeneration, gpuValues.devAbsTol, numThreads, gpuValues.devCConstant, gpuValues.devMarsLaunchCon, gpuValues.devTime_steps, gpuValues.devY_steps, gpuValues.devGamma_steps, gpuValues.devTau_steps, gpuValues.devAccel_steps, gpuValues.devFuel_steps);
        // cudaEventRecord(kernelEnd);

        rk4CPUSim(individual, absTolInput, cConstant, marsLaunchCon, time_steps, y_steps, gamma_steps, tau_steps, accel_steps, fuel_steps);

        // copy the result of the kernel onto the host
        // cudaMemcpy(generation, gpuValues.devGeneration, numThreads * sizeof(Child), cudaMemcpyDeviceToHost);

        //Output how many simulation cycles the code has done
        //std::cout << "\nCompleted simulation cycle " << simulationNum << " for this generation.\n";
        std::cout << "_";

        //simulationNum++; //Add one to simulation num for the next cycle

        // float kernelT;
        
        // cudaEventSynchronize(kernelEnd);

        // cudaEventElapsedTime(&kernelT, kernelStart, kernelEnd);
        
        // calcPerS = numThreads / (kernelT / 1000.0); // how many times the Runge Kutta algorithm ran in the kernel per second


        //Go through simulated individuals to see how many are finished
        // for (int i = 0; i < cConstant->num_individuals; i++)
        // {
        //     if(generation[i].simStatus == COMPLETED_SIM) {
        //         completedNum++;
        //     }
        // }
    } while (individual.simStatus != COMPLETED_SIM);
} */

// Called by optimize() in optimization.cu
//void callRK(const int numThreads, const int blockThreads, Child *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, const cudaConstants* cConstant, PlanetInfo *marsLaunchCon) {
void callRK(Child & individual, double absTolInput, const cudaConstants* cConstant, elements<double> *marsLaunchCon, std::vector<double> & time_steps, std::vector<elements<double>> & y_steps, std::vector<double> & gamma_steps, std::vector<double> & tau_steps, std::vector<double> & accel_steps, std::vector<double> & fuel_steps, double timeInitial) {

    //which (if any) of these do we need to reset here?
    individual.simStatus = INITIAL_SIM;
    individual.errorStatus = NOT_RUN;
    individual.simStartTime = 0;
    individual.stepCount = 0;
    
    double stepSize = (cConstant->minSimVals[TRIPTIME_OFFSET] - timeInitial)/cConstant->max_numsteps;
    
    //Do while loop will keep simulating while there are still children to simulate
    do {

        rk4CPUSim(individual, absTolInput, cConstant, marsLaunchCon, time_steps, y_steps, gamma_steps, tau_steps, accel_steps, fuel_steps);

        std::cout << "_";
    } while (individual.simStatus != COMPLETED_SIM);

    //Check to see if nans are generated in the finalPos elements
    if (isnan(individual.finalPos.r) ||
            isnan(individual.finalPos.theta) ||
            isnan(individual.finalPos.z) ||
            isnan(individual.finalPos.vr) ||
            isnan(individual.finalPos.vtheta) ||
            isnan(individual.finalPos.vz)  ) {
        //Mark the child with the nan variables with the nan_error flag
        individual.errorStatus = NAN_ERROR;
    }//if it is not a nan, the status has already been made valid or sun_error in rk4CUDASim

    //Now that the status has been determined, there is enough information to set pos and speed diffs
    //The two functions will look at the child's errorStatus and set the diffs based on that
    individual.getPosDiff(cConstant);
    individual.getSpeedDiff(cConstant);
    individual.getHorzVelDiff(cConstant);
    individual.getVertVelDiff(cConstant);
    individual.getOrbitPosDiff(cConstant);
    individual.getOrbitSpeedDiff(cConstant);

    for (int j = 0; j < cConstant->missionObjectives.size(); j++) {
        individual.objTargetDiffs[j] = abs(individual.getParameters(cConstant->missionObjectives[j]) - cConstant->missionObjectives[j].target);
    }

    //Get progress now that all of the outputs are ready
    individual.getProgress(cConstant);
        
}

// Called by optimize() in optimization.cu
// This should be used when we don't need a MATLAB output
void callRKBasic(Child& individual, double absTolInput, const cudaConstants* cConstant, elements<double> *marsLaunchCon, double timeInitial) {
    std::vector<double> time_steps(0,0);
    std::vector<elements<double>> y_steps(0,0);
    std::vector<double> gamma_steps(0,0); 
    std::vector<double> tau_steps(0,0); 
    std::vector<double> accel_steps(0,0);
    std::vector<double> fuel_steps(0,0);
    callRK(individual, absTolInput, cConstant, marsLaunchCon, time_steps, y_steps, gamma_steps, tau_steps, accel_steps, fuel_steps, timeInitial); //this doesn't work? what should i be doing
        
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
void rk4CPUSim(Child& individual, double absTolInput, const cudaConstants* cConstant, elements<double> *marsLaunchCon, std::vector<double>& time_steps, std::vector<elements<double>>& y_steps, std::vector<double>& gamma_steps, std::vector<double>& tau_steps, std::vector<double>& accel_steps, std::vector<double>& fuel_steps) {
   
    // individual.simNum++;
    //Check if this child has been simulated already
    if (individual.simStatus != COMPLETED_SIM) {
        //If not, run it's simulation

        //Check to see if the child is about to be simulated too many times
        if ((individual.simNum + 1) > cConstant->maxSimNum) {

            //Assign an error to the child because it has been running for too long
            individual.errorStatus = SIMNUM_ERROR;

            //Mark it as having completed it's run
            individual.simStatus = COMPLETED_SIM;

            //Quit the simulation
            return;
        }

        //If not, run the next simulation and increase the sim num to reflect this
        individual.simNum++;
        
        rkParameters<double> rkParameters = individual.startParams; // get the parameters for this thread

        // storing copies of the input values
        double stepSize;
        double startTime;
        double endTime;
        double curAccel = 0;

        //Stores the total mass of the fuel spent during the simulation
        double massFuelSpent;

        //Initial time and position of the child
        double curTime;
        elements<double> curPos;

        //Stores the angular momentum if the individual is entering an SOI
        //  Used to check if the assist was bad
        double soiEntryh;

        //Set the initial curTime and curPos depending on if the child has been ran
        if (individual.simStatus == INITIAL_SIM) {
            //Child has not been simulated, set the initial curTime to the start time of the simulation
            //Set at 0 initially as no time has passed within the simulation yet
            curTime = 0;
            //Set the start time to the total trip time
            startTime = curTime;
            // start with the initial conditions of the spacecraft
            curPos = rkParameters.y0; 

            massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0
        }
        else {
            //Child has been partially simulated (means it has entered a SOI), set initial curTime to the child's simStartTime variable
            //That will be set to the simStartTime variable for child
            curTime = individual.simStartTime;
            //Set the start time to this simulation's start time
            startTime = individual.simStartTime;
            //Get the child's simStartPos, which will have the elements of the child at the last step of the last simulation
            curPos = individual.simStartPos;

            //Get the mass of how much fuel the child has spent on previous simulations
            massFuelSpent = individual.fuelSpent;
        }

        //Check to see if this simulation occurs within a sphere of influence
        //This determines the endTime
        if(individual.simStatus == INSIDE_SOI) {
            //If this is a mission inside a SOI, set the final to be a gravAssistTime difference from the start time
            endTime = startTime + cConstant->gravAssistTime;
            //endTime = startTime + ((rkParameters.tripTime - startTime) * cConstant->gravAssistTimeFrac);

            //Calculate the entry angular momentum
            soiEntryh = rkParameters.y0.r * rkParameters.y0.vtheta;

            //Check to make sure that endTime is not further than tripTime
            if(endTime > rkParameters.tripTime) {
                endTime = rkParameters.tripTime;
            }
        }
        //In all other scenerios, the end time is the triptime
        else {
            endTime = rkParameters.tripTime;
        }

        thruster<double> thrust(cConstant);

        bool coast; //=1 means thrusting, from calc_coast()

        elements<double> error; // holds output of previous value from rkCalc

        //Calculation of the simulation's step size
        //  Setting the value to be as small as possible
        stepSize = (endTime - startTime) / cConstant->max_numsteps;

        while (curTime < endTime) {
            //Save the time of the current step
            time_steps[individual.stepCount] = curTime;

            // Check the thruster type before performing calculations
            if (cConstant->thruster_type == thruster<double>::NO_THRUST) {
                coast = curAccel = 0;
            }
            else {
                //Calc if the thruster is activated
                coast = calc_coast(rkParameters.coeff, curTime, rkParameters.tripTime, thrust);
                
                //Calc current step's acceleration
                curAccel = calc_accel(curPos.r, curPos.z, thrust, massFuelSpent, stepSize, coast, static_cast<double>(cConstant->wet_mass), cConstant);
                
                //Record the gamma value of the current step
                gamma_steps.push_back(calc_gamma(individual.startParams.coeff, curTime, rkParameters.tripTime));

                //Record the tau value of the current step
                tau_steps.push_back(calc_tau(individual.startParams.coeff, curTime, rkParameters.tripTime));

                //Update the child with how much fuel it has used 
                individual.fuelSpent = massFuelSpent;
            }
            //Record the acceleration and fuel spent values of the current step
            accel_steps[individual.stepCount] = curAccel;
            fuel_steps[individual.stepCount] = massFuelSpent;

            //Needs to be triptime - curtime to get the correct index for mars
            //when curtime = triptime, this will give us the final position of mars at impact
            //this is because getConditionDev takes in seconds before the spacecraft reaches the target
            elements<double> mars = getConditionDev(rkParameters.tripTime - curTime, cConstant, marsLaunchCon);

            //calculate the distance between mars and the spacecraft (|R|^2)
            double marsCraftDist = sqrt(pow(mars.r, 2) + pow(curPos.r, 2) + pow(curPos.z - mars.z, 2) - (2*curPos.r*mars.r*cos(mars.theta-curPos.theta)));

            //See if the child is closest it has been to Mars so far this run
            //This is only updated if Mars is in betweem the craft and the target
            if (marsCraftDist < individual.minMarsDist) {
                individual.minMarsDist = marsCraftDist;
            }
            
            // calculate k values and get new value of y
            rkCalc(curTime, rkParameters.tripTime, stepSize, curPos, rkParameters.coeff, curAccel, error, mars, marsCraftDist); 

            curTime += stepSize; // update the current time in the simulation
            
            //TODO: This is unused as of summer 2023 (as the minimum step size was nearly always exceeded), should this be a permenant change?
            //stepSize *= calc_scalingFactor(curPos-error,error,absTol, cConstant->doublePrecThresh); // Alter the step size for the next iteration

            //// The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            //if (stepSize > (endTime - startTime) / cConstant->min_numsteps) {
            //    stepSize = (endTime - startTime) / cConstant->min_numsteps;
            //}
            //else if (stepSize < (endTime - startTime) / cConstant->max_numsteps){
            //    stepSize = (endTime - startTime) / cConstant->max_numsteps;
            //}

            if ( (curTime + stepSize) > endTime) {
                stepSize = (endTime - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if ( sqrt(pow(curPos.r,2) + pow(curPos.z,2)) < cConstant->sun_r_min) { //maybe issue is with using pow? I doubt it, but we could always try curPos.r*curPos.r + curPos.z*curPos.z < sun_r_min*sun_r_min?
                //This is a bad result, needs to be set to be removed
                //Setting the child's status to be a sun error
                individual.errorStatus = SUN_ERROR;//Are all the children's errorStatus set to SUN_ERROR?

                //Set the child's sim status to complete so it ins't ran any further
                individual.simStatus = COMPLETED_SIM; 

                return;
            }
            
            //Check to see if the child is too close to Mars
            if(marsCraftDist < cConstant->gravAssistDist){
                //This is a bad result, needs to be set to be removed
                //Setting the child's status to be a mars error
                individual.errorStatus = MARS_ERROR;
                
                //Set the child's sim status to complete so it ins't ran any further
                individual.simStatus = COMPLETED_SIM;

                return;
            }

            //Check to see if the curTime is less than triptime
            if (curTime < endTime) {
                //If so, check to see if the child triggers the conditions for a new run

                //See if child has entered MSOI
                if (marsCraftDist < MSOI*cConstant->MSOI_scale && individual.simStatus != INSIDE_SOI) {
                    //Set the final time and position for this portion of the child's simulation
                    individual.simStartTime = curTime;
                    individual.simStartPos = curPos;

                    //Set the child status to inside SOI
                    individual.simStatus = INSIDE_SOI;

                    return;
                }

                //Check if child has exited MSOI after being inside it
                if (marsCraftDist > MSOI*cConstant->MSOI_scale && individual.simStatus == INSIDE_SOI) {
                    //Set the final time and position for this portion of the child's simulation
                    individual.simStartTime = curTime;
                    individual.simStartPos = curPos;

                    //Set the child status to outide SOI
                    individual.simStatus = OUTSIDE_SOI;

                    //Calculate the change in angular momentum during the child's time in the SOI
                    individual.orbithChange = (curPos.r * curPos.vtheta) - soiEntryh;

                    //Check to make sure the orbithChange isn't less than 0, if the value is less than 0 the intercept calculations will mess up
                    if (individual.orbithChange < 1e-14) {
                        individual.orbithChange = 5e-15;
                    }

                    //Check to see if this was a bad assist
                    // if (curPos.r * curPos.vtheta < soiEntryh) {
                    //     //If it is a bad assist, set the error
                    //     individual.errorStatus = BAD_ASSIST;
                    //     individual.simStatus = COMPLETED_SIM;
                    // }

                    return;
                }

                //Save the position and velocity of the current time step
                y_steps[individual.stepCount] = curPos;
            
                //count the steps taken for this threads calculations
                individual.stepCount++;

            }  //end if (curTime < endTime)
        } // end while (curTime < endTime)

        //Check to see if this simulation has completed its total runtime
        if (endTime >= rkParameters.tripTime) {
            //If here, the individual's simulation has completely ended, regardless of its current simStatus
            
            //Give the child its final calculated position
            individual.finalPos = curPos;

            //if it is not a SUN_ERROR then it is valid
            individual.errorStatus = VALID;

            //The child has finished it's run, set the simStatus to completed
            individual.simStatus = COMPLETED_SIM;

            //individual.fuelSpent = massFuelSpent;
        }
        //else: endTime != tripTime
        else {
            //This means the simulation was for an individual in a SOI, but it never escaped the SOI in the estimated time
            //The individual will need be run through a SOI simulation again 
            //The sim status would already be INSIDE_SOI, so only the position and time need to be recorded
            individual.simStartPos = curPos;
            individual.simStartTime = curTime;  
        }

        return;
        
    }
}

template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
                                   T stepSize, elements<T> & y_new, const T & absTol, const cudaConstants * cConstants) {
    // Set the first element of the solution vector to the conditions of the Planet on impact date (Oct. 5, 2022 for Earth)
    y_new = y0;
    elements<T> error;
    T curTime = timeFinal; // setting time equal to the start time

    while( curTime > timeInitial) {  // iterates in reverse
        //calculate k values
        rkCalcPlanet(curTime, timeFinal, stepSize, y_new, error);

        //array of time output as t         
        curTime -= stepSize;

        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(y_new-error,error,absTol, cConstants->doublePrecThresh);

        //The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize > (timeFinal-timeInitial) / cConstants->min_numsteps) {
            stepSize = (timeFinal-timeInitial) / cConstants->min_numsteps;
        }
        else if (stepSize < ((timeFinal-timeInitial) / cConstants->max_numsteps)) {
            stepSize = (timeFinal-timeInitial) / cConstants->max_numsteps;
        }

        // shorten the last step to end exactly at time final
        if ( (curTime-stepSize) < timeInitial) {
            stepSize = curTime-timeInitial;
        }
    } //end of while
}

template <class T>  void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, coefficients<T> & coeff, const T & accel, 
                                                    elements<T> & error, const elements<T> & mars, T & marsCraftDist) {

    // k variables for Runge-Kutta calculation of y_new
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    // Coefficients from MATLAB's implementation of ode45
    // Our calculation of k has the time step built into it (see motion_equations.cpp)
    //marsCraftDist = calcMarsCraftDist(y_new, mars);
    k1 = calc_k(stepSize, y_new, coeff, accel, curTime, timeFinal, mars, marsCraftDist); 
    k2 = calc_k(stepSize, y_new+k1*(static_cast <double> (1)/static_cast <double> (5)), coeff, accel, curTime+((static_cast <double> (1)/static_cast <double> (5))*stepSize), timeFinal, mars, marsCraftDist); 
    k3 = calc_k(stepSize, y_new+k1*(static_cast <double> (3)/static_cast <double> (40))+k2*(static_cast <double> (9)/static_cast <double> (40)), coeff, accel, curTime+((static_cast <double> (3)/static_cast <double> (10))*stepSize), timeFinal, mars, marsCraftDist);   
    k4 = calc_k(stepSize, y_new+k1*(static_cast <double> (44)/static_cast <double> (45))+k2*(static_cast <double> (-56)/static_cast <double> (15))+k3*(static_cast <double> (32)/static_cast <double> (9)), coeff, accel, curTime+((static_cast <double> (4)/static_cast <double> (5))*stepSize), timeFinal, mars, marsCraftDist); 
    k5 = calc_k(stepSize, y_new+k1*(static_cast <double> (19372)/static_cast <double> (6561))+k2*(static_cast <double> (-25360)/static_cast <double> (2187))+k3*(static_cast <double> (64448)/static_cast <double> (6561))+k4*(static_cast <double> (-212)/static_cast <double> (729)), coeff, accel, curTime+((static_cast <double> (8)/static_cast <double> (9))*stepSize), timeFinal, mars, marsCraftDist); 
    k6 = calc_k(stepSize, y_new+k1*(static_cast <double> (9017)/static_cast <double> (3168))+k2*(static_cast <double> (-355)/static_cast <double> (33))+k3*(static_cast <double> (46732)/static_cast <double> (5247))+k4*(static_cast <double> (49)/static_cast <double> (176))+k5*(static_cast <double> (-5103)/static_cast <double> (18656)), coeff, accel, curTime+stepSize, timeFinal, mars, marsCraftDist);  
    k7 = calc_k(stepSize, y_new+k1*(static_cast <double> (35)/static_cast <double> (384))+k3*(static_cast <double> (500)/static_cast <double> (1113))+k4*(static_cast <double> (125)/static_cast <double> (192))+k5*(static_cast <double> (-2187)/static_cast <double> (6784))+k6*(static_cast <double> (11)/static_cast <double> (84)), coeff, accel, curTime+stepSize, timeFinal, mars, marsCraftDist);  

    // New value
    y_new = y_new + k1*(static_cast <double> (35)/static_cast <double> (384)) + k3*(static_cast <double> (500)/static_cast <double> (1113)) + k4*(static_cast <double> (125)/static_cast <double> (192)) - k5*(static_cast <double> (2187)/static_cast <double> (6784)) + k6*(static_cast <double> (11)/static_cast <double> (84)) + k7*(static_cast <double> (0)/static_cast <double> (40));  

    // Error 
    // See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    // Dormand-Prince : no error between GPU and CPU
    //y_prev = k1*5179./57600 + k3*7571./16695 + k4*393./640 - k5*92097./339200 + k6*187./2100 + k7*1./40;  
    //error = y_new-y_prev;

    // MATLAB code : ERROR between GPU and CPU
    //error = k1*(71)/(57600) + k3*(-71)/(16695) + k4*(71)/(1920)
    //- k5*(17253)/(339200) + k6*(22)/(525) + k7*(-1)/(40);

    // (outdated comment, but still important) Without k7 : no error between GPU and CPU (this has to do with the comment below)
    
    // Comonents of error are going to be really small. Need to make sure they are not too small to do anything with in calc_scalingFactor
    error =  ((k1*(static_cast <double> (71)/static_cast <double> (57600))) + (k3*(static_cast <double> (-71)/static_cast <double> (16695))) + (k4*(static_cast <double> (71)/static_cast <double> (1920)))  + (k5*(static_cast <double> (-17253)/static_cast <double> (339200))) + (k6*(static_cast <double> (22)/static_cast <double> (525)))) + (k7*(static_cast <double> (-1)/static_cast <double> (40)));
}

// The stepSize value that is inputted is assumed to be a positive value
template <class T> void rkCalcPlanet(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, elements<T> & error) {
    // Runge-Kutta algorithm    
    // k variables for Runge-Kutta calculation of y_new
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    
    stepSize *= -1; // Make this copy of stepSize negative as it goes backwards

    //calc_k multiplies all values by the stepSize internally.
    k1 = calc_kPlanet(stepSize, y_new, curTime, timeFinal);        
    k2 = calc_kPlanet(stepSize, y_new+(k1*(static_cast <double> (1)/static_cast <double> (5))), curTime+((static_cast <double> (1)/static_cast <double> (5))*stepSize), timeFinal);   
    k3 = calc_kPlanet(stepSize, y_new+(k1*(static_cast <double> (3)/static_cast <double> (40)))+(k2*(static_cast <double> (9)/static_cast <double> (40))), curTime+((static_cast <double> (3)/static_cast <double> (10))*stepSize), timeFinal);   
    k4 = calc_kPlanet(stepSize, y_new+(k1*(static_cast <double> (44)/static_cast <double> (45)))+(k2*(static_cast <double> (-56)/static_cast <double> (15)))+(k3*(static_cast <double> (32)/static_cast <double> (9))), curTime+((static_cast <double> (4)/static_cast <double> (5))*stepSize), timeFinal);    
    k5 = calc_kPlanet(stepSize, y_new+(k1*(static_cast <double> (19372)/static_cast <double> (6561)))+(k2*(static_cast <double> (-25360)/static_cast <double> (2187)))+(k3*(static_cast <double> (64448)/static_cast <double> (6561)))+(k4*(static_cast <double> (-212)/static_cast <double> (729))), curTime+((static_cast <double> (8)/static_cast <double> (9))*stepSize), timeFinal);        
    k6 = calc_kPlanet(stepSize, y_new+(k1*(static_cast <double> (9017)/static_cast <double> (3168)))+(k2*(static_cast <double> (-355)/static_cast <double> (33)))+(k3*(static_cast <double> (46732)/static_cast <double> (5247)))+(k4*(static_cast <double> (49)/static_cast <double> (176)))+(k5*(static_cast <double> (-5103)/static_cast <double> (18656))), curTime+stepSize, timeFinal);        
    k7 = calc_kPlanet(stepSize, y_new+(k1*(static_cast <double> (35)/static_cast <double> (384)))+(k3*(static_cast <double> (500)/static_cast <double> (1113)))+(k4*(static_cast <double> (125)/static_cast <double> (192)))+(k5*(static_cast <double> (-2187)/static_cast <double> (6784)))+(k6*(static_cast <double> (11)/static_cast <double> (84))), curTime+stepSize, timeFinal);  

    //Error 
    //See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    //v = y_new + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;  

    //New value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    y_new = y_new + (k1*(static_cast <double> (35)/static_cast <double> (384))) + (k3*(static_cast <double> (500)/static_cast <double> (1113))) + (k4*(static_cast <double> (125)/static_cast <double> (192))) + (k5*(static_cast <double> (-2187)/static_cast <double> (6784))) + (k6*(static_cast <double> (11)/static_cast <double> (84)));  

    // Error 
    // See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    // Dormand-Prince : no error between GPU and CPU
    //y_prev = k1*5179./57600 + k3*7571./16695 + k4*393./640 - k5*92097./339200 + k6*187./2100 + k7*1./40;  
    //error = y_new-y_prev;

    error = (k1*(static_cast <double> (71)/static_cast <double> (57600))) + (k3*(static_cast <double> (-71)/static_cast <double> (16695))) + (k4*(static_cast <double> (71)/static_cast <double> (1920))) - (k5*(static_cast <double> (17253)/static_cast <double> (339200))) + (k6*(static_cast <double> (22)/static_cast <double> (525))) + (k7*(static_cast <double> (-1)/static_cast <double> (40)));    
}

template <class T>  T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, const double precThresh) {
    // relative total error is the total error of all coponents of y which is used in scale.
    // scale is used to determine the next step size.
    T normTotError, scale;

    // relative error (unitless) 
    elements<T> pmError(difference.r/previous.r, difference.theta/previous.theta, difference.z/previous.z, 
    difference.vr/previous.vr,  difference.vtheta/previous.vtheta, difference.vz/previous.vz);

    if (!pmLimitCheck(pmError, precThresh)) {
        // pmError is too small!
        // Keep time step the same
        // Complicated rational:
        // If error is so small the error is at the double error limit
        //     - Could mean time step needs to be smaller, thus more presice
        //     - Or could mean that time step needs to be longer to get a reasonable error
        return 1.0; // Change step size by 0%
    }

    // elements<T> pmError(previous.r, previous.theta, previous.z, previous.vr,  previous.vtheta, previous.vz);

    // square root of sum of squares of the error from the 6 elements to determine the scale for the time step of the next iteration
    normTotError = sqrt(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2));
    scale = pow((absTol/normTotError),0.2);

    return scale;   
}

template <class T>  bool pmLimitCheck(const elements<T> & pmError, const double precThresh){
    //It is possible this is a major resource drain. This might be faster to square everything and not use fabs (floating point abs)
    if( (fabs(pmError.r) < precThresh ) ||
        (fabs(pmError.theta) < precThresh ) ||
        (fabs(pmError.z) < precThresh ) ||
        (fabs(pmError.vr) < precThresh ) ||
        (fabs(pmError.vtheta) < precThresh ) ||
        (fabs(pmError.vz) < precThresh ) )
    {
        //Error it too small for precise calculation of step size
        return false;
    } else {
        //Error is large enough that accurate step size can be computed
        // All error values are at least within 14 orders of magnitued of their original guesses
        return true;
    }
}