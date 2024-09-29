#include "mex.h"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <fstream>  // Include fstream for file I/O

// Define data structures
struct DendriteParams {
    double alpha;
    double b;
    double Tau;
    double TauR;
    double NaX;
    double gc;
    int ID;
};

struct StateVariables {
    double u;
    double v;
};

struct CouplingEntry {
    int row;
    int col;
    double value;
};

// Function prototypes
void RK4Solver(
    const std::vector<DendriteParams>& dendrites,
    const std::vector<CouplingEntry>& couplingData,
    const std::unordered_map<int, int>& dendriteIDtoIndex,
    const std::unordered_map<int, int>& stimDendriteIDtoIndex,
    const std::vector<std::vector<double>>& StimuliMatrix,
    const std::vector<int>& stimDendriteIDs,
    std::vector<std::vector<double>>& solution,
    double h,
    int numTimeSteps,
    const std::vector<double>& t,
    const std::vector<StateVariables>& initialStates
);

void computeRHS(
    const std::vector<StateVariables>& states,
    const std::vector<DendriteParams>& dendrites,
    const std::vector<CouplingEntry>& couplingData,
    const std::unordered_map<int, int>& dendriteIDtoIndex,
    std::vector<double>& couplingTerm,
    std::vector<double>& Stim,
    std::vector<StateVariables>& derivatives,
    const std::unordered_map<int, int>& stimDendriteIDtoIndex,
    const std::vector<std::vector<double>>& StimuliMatrix,
    int timeIndex
);

// Entry point
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check the number of inputs and outputs
    if (nrhs != 8) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidNumInputs", "Eight inputs required.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidNumOutputs", "Two outputs required.");
    }

    // Parse inputs
    // 0: paramsArray (numDendrites x 7)
    // 1: couplingData (numCouplings x 3)
    // 2: StimuliMatrix (numStimDendrites x numTimeSteps)
    // 3: stimDendriteIDArray (numStimDendrites x 1)
    // 4: X0 (numStateVars x 1)
    // 5: h (scalar)
    // 6: numTimeSteps (scalar)
    // 7: dendriteIDArray (numDendrites x 1)

    // 1. Parse paramsArray
    const mxArray* paramsArray_mx = prhs[0];
    mwSize numDendrites = mxGetM(paramsArray_mx);
    double* paramsArray_ptr = mxGetPr(paramsArray_mx);
    if (paramsArray_ptr == nullptr) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "paramsArray must be of type double.");
    }
    std::vector<DendriteParams> dendrites(numDendrites);

    for (mwSize i = 0; i < numDendrites; ++i) {
        dendrites[i].alpha = paramsArray_ptr[i];
        dendrites[i].b     = paramsArray_ptr[i + numDendrites];
        dendrites[i].Tau   = paramsArray_ptr[i + 2 * numDendrites];
        dendrites[i].TauR  = paramsArray_ptr[i + 3 * numDendrites];
        dendrites[i].NaX   = paramsArray_ptr[i + 4 * numDendrites];
        dendrites[i].gc    = paramsArray_ptr[i + 5 * numDendrites];
        dendrites[i].ID    = static_cast<int>(paramsArray_ptr[i + 6 * numDendrites]);
    }

    // 2. Parse couplingData
    const mxArray* couplingData_mx = prhs[1];
    mwSize numCouplings = mxGetM(couplingData_mx);
    double* couplingData_ptr = mxGetPr(couplingData_mx);
    if (couplingData_ptr == nullptr) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "couplingData must be of type double.");
    }
    std::vector<CouplingEntry> couplingData(numCouplings);

    for (mwSize i = 0; i < numCouplings; ++i) {
        couplingData[i].row   = static_cast<int>(couplingData_ptr[i]);                         // row index (zero-based)
        couplingData[i].col   = static_cast<int>(couplingData_ptr[i + numCouplings]);          // col index (zero-based)
        couplingData[i].value = couplingData_ptr[i + 2 * numCouplings];                        // value
    }

    // 3. Parse StimuliMatrix
    const mxArray* StimuliMatrix_mx = prhs[2];
    mwSize numStimDendrites = mxGetM(StimuliMatrix_mx);
    mwSize numTimeSteps = mxGetN(StimuliMatrix_mx);
    double* StimuliMatrix_ptr = mxGetPr(StimuliMatrix_mx);
    if (StimuliMatrix_ptr == nullptr) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "StimuliMatrix must be of type double.");
    }

    std::vector<std::vector<double>> StimuliMatrix(numStimDendrites, std::vector<double>(numTimeSteps));
    for (mwSize i = 0; i < numStimDendrites; ++i) {
        for (mwSize j = 0; j < numTimeSteps; ++j) {
            StimuliMatrix[i][j] = StimuliMatrix_ptr[i + j * numStimDendrites];
        }
    }

    // 4. Parse stimDendriteIDArray
    const mxArray* stimDendriteIDArray_mx = prhs[3];
    double* stimDendriteIDArray_ptr = mxGetPr(stimDendriteIDArray_mx);
    if (stimDendriteIDArray_ptr == nullptr) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "stimDendriteIDArray must be of type double.");
    }
    std::vector<int> stimDendriteIDs(numStimDendrites);
    for (mwSize i = 0; i < numStimDendrites; ++i) {
        stimDendriteIDs[i] = static_cast<int>(stimDendriteIDArray_ptr[i]);
    }

    // 5. Parse X0
    const mxArray* X0_mx = prhs[4];
    double* X0_ptr = mxGetPr(X0_mx);
    if (X0_ptr == nullptr) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "X0 must be of type double.");
    }
    mwSize numStateVars = mxGetM(X0_mx); // Should be numDendrites * 2
    if (numStateVars != numDendrites * 2) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "Size of X0 does not match number of dendrites.");
    }
    std::vector<StateVariables> initialStates(numDendrites);
    for (mwSize i = 0; i < numDendrites; ++i) {
        initialStates[i].u = X0_ptr[2 * i];
        initialStates[i].v = X0_ptr[2 * i + 1];
    }

    // 6. Parse h
    double h = mxGetScalar(prhs[5]);

    // 7. Parse numTimeSteps (overrides the one from StimuliMatrix)
    int numTimeSteps_input = static_cast<int>(mxGetScalar(prhs[6]));
    if (numTimeSteps != static_cast<mwSize>(numTimeSteps_input)) {
        mexErrMsgIdAndTxt("mexRK4Solver:timeStepMismatch", "numTimeSteps does not match StimuliMatrix dimensions.");
    }

    // 8. Parse dendriteIDArray
    const mxArray* dendriteIDArray_mx = prhs[7];
    double* dendriteIDArray_ptr = mxGetPr(dendriteIDArray_mx);
    if (dendriteIDArray_ptr == nullptr) {
        mexErrMsgIdAndTxt("mexRK4Solver:invalidInput", "dendriteIDArray must be of type double.");
    }
    std::unordered_map<int, int> dendriteIDtoIndex;
    for (mwSize i = 0; i < numDendrites; ++i) {
        int ID = static_cast<int>(dendriteIDArray_ptr[i]);
        dendriteIDtoIndex[ID] = static_cast<int>(i); // zero-based indexing
    }

    // Create mapping from stimDendriteID to index in StimuliMatrix
    std::unordered_map<int, int> stimDendriteIDtoIndex;
    for (mwSize i = 0; i < numStimDendrites; ++i) {
        stimDendriteIDtoIndex[stimDendriteIDs[i]] = static_cast<int>(i); // zero-based indexing
    }

    // ** Start of Logging Code **

    // Open a log file to write the data
    std::ofstream logFile("mexRK4Solver_log.txt");
    if (!logFile.is_open()) {
        mexErrMsgIdAndTxt("mexRK4Solver:logFileError", "Unable to open log file for writing.");
    }

    // Write data to the log file

    // 1. Write dendrite parameters
    logFile << "Dendrite Parameters (paramsArray):\n";
    logFile << "Index\talpha\tb\tTau\tTauR\tNaX\tgc\tID\n";
    for (mwSize i = 0; i < numDendrites; ++i) {
        logFile << i << "\t"
                << dendrites[i].alpha << "\t"
                << dendrites[i].b << "\t"
                << dendrites[i].Tau << "\t"
                << dendrites[i].TauR << "\t"
                << dendrites[i].NaX << "\t"
                << dendrites[i].gc << "\t"
                << dendrites[i].ID << "\n";
    }
    logFile << "\n";

    // 2. Write coupling data
    logFile << "Coupling Data (couplingData):\n";
    logFile << "Entry\tRow\tCol\tValue\n";
    for (mwSize i = 0; i < numCouplings; ++i) {
        logFile << i << "\t"
                << couplingData[i].row << "\t"
                << couplingData[i].col << "\t"
                << couplingData[i].value << "\n";
    }
    logFile << "\n";

    // 3. Write Stimuli Matrix
    logFile << "Stimuli Matrix (StimuliMatrix):\n";
    logFile << "Stimulus Index (Row in StimuliMatrix) corresponds to Dendrite ID in stimDendriteIDs\n";
    for (mwSize i = 0; i < numStimDendrites; ++i) {
        logFile << "Stimulus Index: " << i << ", Dendrite ID: " << stimDendriteIDs[i] << "\n";
        logFile << "Time Step\tStimulus Value\n";
        for (mwSize j = 0; j < numTimeSteps; ++j) {
            logFile << j << "\t" << StimuliMatrix[i][j] << "\n";
        }
        logFile << "\n";
    }
    logFile << "\n";

    // 4. Write stimDendriteIDtoIndex mapping
    logFile << "Stimulated Dendrite ID to Index Mapping (stimDendriteIDtoIndex):\n";
    logFile << "Dendrite ID\tStimulus Index\n";
    for (auto it = stimDendriteIDtoIndex.begin(); it != stimDendriteIDtoIndex.end(); ++it) {
        logFile << it->first << "\t" << it->second << "\n";
    }
    logFile << "\n";

    // 5. Write dendriteIDtoIndex mapping
    logFile << "Dendrite ID to Index Mapping (dendriteIDtoIndex):\n";
    logFile << "Dendrite ID\tIndex\n";
    for (auto it = dendriteIDtoIndex.begin(); it != dendriteIDtoIndex.end(); ++it) {
        logFile << it->first << "\t" << it->second << "\n";
    }
    logFile << "\n";

    // 6. Write initial conditions
    logFile << "Initial Conditions (X0):\n";
    logFile << "Index\tu\tv\n";
    for (mwSize i = 0; i < numDendrites; ++i) {
        logFile << i << "\t" << initialStates[i].u << "\t" << initialStates[i].v << "\n";
    }
    logFile << "\n";

    // Close the log file
    logFile.close();

    // ** End of Logging Code **

    // Prepare time vector t
    std::vector<double> t(numTimeSteps);
    for (int i = 0; i < numTimeSteps; ++i) {
        t[i] = i * h;
    }

    // Prepare solution matrix
    std::vector<std::vector<double>> solution(numStateVars, std::vector<double>(numTimeSteps));

    // Run the RK4 solver
    RK4Solver(dendrites, couplingData, dendriteIDtoIndex, stimDendriteIDtoIndex, StimuliMatrix, stimDendriteIDs, solution, h, numTimeSteps, t, initialStates);

    // Prepare outputs
    // Output t (time vector)
    plhs[0] = mxCreateDoubleMatrix(1, numTimeSteps, mxREAL);
    double* t_out = mxGetPr(plhs[0]);
    std::copy(t.begin(), t.end(), t_out);

    // Output solution matrix
    plhs[1] = mxCreateDoubleMatrix(numStateVars, numTimeSteps, mxREAL);
    double* solution_out = mxGetPr(plhs[1]);
    for (mwSize i = 0; i < numStateVars; ++i) {
        for (int j = 0; j < numTimeSteps; ++j) {
            solution_out[i + j * numStateVars] = solution[i][j];
        }
    }
}

// Implement the RK4 solver
void RK4Solver(
    const std::vector<DendriteParams>& dendrites,
    const std::vector<CouplingEntry>& couplingData,
    const std::unordered_map<int, int>& dendriteIDtoIndex,
    const std::unordered_map<int, int>& stimDendriteIDtoIndex,
    const std::vector<std::vector<double>>& StimuliMatrix,
    const std::vector<int>& stimDendriteIDs,
    std::vector<std::vector<double>>& solution,
    double h,
    int numTimeSteps,
    const std::vector<double>& t,
    const std::vector<StateVariables>& initialStates
)
{
    size_t numDendrites = dendrites.size();
    // Initialize state variables
    std::vector<StateVariables> states = initialStates;
    // Initialize derivatives
    std::vector<StateVariables> k1(numDendrites), k2(numDendrites), k3(numDendrites), k4(numDendrites);
    std::vector<StateVariables> tempStates(numDendrites);

    // Initialize couplingTerm and Stim
    std::vector<double> couplingTerm(numDendrites, 0.0);
    std::vector<double> Stim(numDendrites, 0.0);

    // Store initial conditions in the solution
    for (mwSize i = 0; i < numDendrites; ++i) {
        solution[2 * i][0] = states[i].u;
        solution[2 * i + 1][0] = states[i].v;
    }

    // Time integration loop
    for (int step = 1; step < numTimeSteps; ++step) {
        double currentTime = t[step - 1];

        // Compute k1
        computeRHS(states, dendrites, couplingData, dendriteIDtoIndex, couplingTerm, Stim, k1, stimDendriteIDtoIndex, StimuliMatrix, step - 1);

        // Compute tempStates for k2
        for (mwSize i = 0; i < numDendrites; ++i) {
            tempStates[i].u = states[i].u + 0.5 * h * k1[i].u;
            tempStates[i].v = states[i].v + 0.5 * h * k1[i].v;
        }
        computeRHS(tempStates, dendrites, couplingData, dendriteIDtoIndex, couplingTerm, Stim, k2, stimDendriteIDtoIndex, StimuliMatrix, step - 1);

        // Compute tempStates for k3
        for (mwSize i = 0; i < numDendrites; ++i) {
            tempStates[i].u = states[i].u + 0.5 * h * k2[i].u;
            tempStates[i].v = states[i].v + 0.5 * h * k2[i].v;
        }
        computeRHS(tempStates, dendrites, couplingData, dendriteIDtoIndex, couplingTerm, Stim, k3, stimDendriteIDtoIndex, StimuliMatrix, step - 1);

        // Compute tempStates for k4
        for (mwSize i = 0; i < numDendrites; ++i) {
            tempStates[i].u = states[i].u + h * k3[i].u;
            tempStates[i].v = states[i].v + h * k3[i].v;
        }
        computeRHS(tempStates, dendrites, couplingData, dendriteIDtoIndex, couplingTerm, Stim, k4, stimDendriteIDtoIndex, StimuliMatrix, step);

        // Update states
        for (mwSize i = 0; i < numDendrites; ++i) {
            states[i].u += (h / 6.0) * (k1[i].u + 2 * k2[i].u + 2 * k3[i].u + k4[i].u);
            states[i].v += (h / 6.0) * (k1[i].v + 2 * k2[i].v + 2 * k3[i].v + k4[i].v);

            // Store the new states
            solution[2 * i][step] = states[i].u;
            solution[2 * i + 1][step] = states[i].v;
        }
    }
}

// Implement the RHS computation
void computeRHS(
    const std::vector<StateVariables>& states,
    const std::vector<DendriteParams>& dendrites,
    const std::vector<CouplingEntry>& couplingData,
    const std::unordered_map<int, int>& dendriteIDtoIndex,
    std::vector<double>& couplingTerm,
    std::vector<double>& Stim,
    std::vector<StateVariables>& derivatives,
    const std::unordered_map<int, int>& stimDendriteIDtoIndex,
    const std::vector<std::vector<double>>& StimuliMatrix,
    int timeIndex
)
{
    size_t numDendrites = dendrites.size();
    // Reset coupling terms
    std::fill(couplingTerm.begin(), couplingTerm.end(), 0.0);

    // Compute coupling terms
    for (const auto& entry : couplingData) {
        int i = entry.row; // source index
        int j = entry.col; // target index
        double gc = entry.value;
        couplingTerm[i] += gc * (states[j].u - states[i].u);
    }

    // Compute stimuli
    std::fill(Stim.begin(), Stim.end(), 0.0);
    for (auto it = stimDendriteIDtoIndex.begin(); it != stimDendriteIDtoIndex.end(); ++it) {
        int ID = it->first;
        int idx = it->second;
        int dendriteIdx = dendriteIDtoIndex.at(ID);
        Stim[dendriteIdx] += StimuliMatrix[idx][timeIndex]; // Accumulate stimuli
    }
   /*for (auto it = stimDendriteIDtoIndex.begin(); it != stimDendriteIDtoIndex.end(); ++it) {
        int ID = it->first;
        int idx = it->second;
        int dendriteIdx = dendriteIDtoIndex.at(ID);
        mexPrintf("id: %f idx: %f \n", ID, idx);
        Stim[dendriteIdx] = StimuliMatrix[idx][timeIndex];
        // mexPrintf("Stim: %f", Stim[dendriteIdx]);
    }*/

    // Compute derivatives
    for (mwSize i = 0; i < numDendrites; ++i) {
        const DendriteParams& params = dendrites[i];
        const StateVariables& state = states[i];
        double u = state.u;
        double v = state.v;

        double du_dt = params.Tau * (params.NaX * (u * (u - 1) * (1 - params.alpha * u) - v) + couplingTerm[i] + Stim[i]);
        double dv_dt = params.TauR * params.b * u;

        derivatives[i].u = du_dt;
        derivatives[i].v = dv_dt;
    }
}
