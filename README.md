# Dendritic Tree Model Simulation

This repository contains a MATLAB implementation of a dendritic tree model for simulating neuronal dynamics. The model incorporates dendrites and a soma, with the ability to simulate electrical activity using the Runge-Kutta 4th order (RK4) numerical integration method. A MEX function written in C++ (`mexRK4Solver.cpp`) accelerates the computation for large-scale simulations.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Compilation](#compilation)
- [Usage](#usage)
  - [Creating Dendrites](#creating-dendrites)
  - [Building the Dendritic Tree Model](#building-the-dendritic-tree-model)
  - [Adding Connections](#adding-connections)
  - [Adding Stimuli](#adding-stimuli)
  - [Running Simulations](#running-simulations)
  - [Visualization](#visualization)
- [Project Structure](#project-structure)
- [Changelog](#changelog)
- [Future Plans](#futute-plans) 

## Overview

The dendritic tree model simulates the electrical activity of neurons by modeling dendrites and their connections to a soma. Each dendrite and the soma are represented as compartments with their own set of differential equations. The model allows for complex connectivity patterns and external stimuli applied to individual dendrites.

The simulation uses:

- **MATLAB** for modeling, visualization, and interfacing.
- **C++ MEX Function** (`mexRK4Solver.cpp`) for efficient numerical integration using the RK4 method.

## Features

- **Modular Design**: Easily add or remove dendrites and connections.
- **Customizable Parameters**: Each dendrite can have unique parameters and initial conditions.
- **External Stimuli**: Apply time-varying stimuli to any dendrite.
- **Efficient Computation**: Accelerated simulations using a C++ MEX function.
- **Visualization Tools**: Visualize the dendritic tree in 3D and the simulation results.

## Prerequisites

- **MATLAB** R2017b or later.
- **C++ Compiler** compatible with MATLAB's MEX (e.g., GCC on Linux, Xcode on macOS, Visual Studio on Windows).
- Basic understanding of neuronal modeling and MATLAB programming.

## Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/ForestSquirrel/dendritic-tree-model.git
   cd dendritic-tree-model
   ```

2. **Add the Project Folder to MATLAB Path**

  Open MATLAB and add the project folder to the MATLAB path:
  ```matlab
  addpath(genpath('/path/to/dendritic-tree-model'));
  ```

## Compilation

### Compile the MEX Function

Before running simulations, you need to compile the `mexRK4Solver.cpp` file:

1. **Set Up a Compatible C++ Compiler**

   Ensure that you have a C++ compiler configured with MATLAB:

   ```matlab
   mex -setup cpp

  Follow the prompts to select a compiler if one isn't already configured.

2. **Compile the MEX Function**

  In MATLAB, navigate to the directory containing mexRK4Solver.cpp and run:
  
  ```matlab
  mex mexRK4Solver.cpp
  ```
  This will generate a MEX file (mexRK4Solver.*) appropriate for your platform.

## Usage
### Creating Dendrites
Create dendrite objects using the Dendrite class. You can specify parameters and initial conditions:

```matlab
% Create a dendrite with ID = 1, custom parameters, and initial conditions
params.alpha = 10;
params.b = 2;
params.Tau = 1700;
params.TauR = 1700;
params.NaX = 0.05;
params.gc = 4;
initialConditions = [0.2; 0.2];

dendrite1 = Dendrite(1, params, initialConditions);
```
### Building the Dendritic Tree Model
Instantiate a DendriteTreeModel object and add dendrites:

```matlab
% Create the model
model = DendriteTreeModel();

% Add dendrites to the model
model = model.addDendrite(dendrite1);
% Add more dendrites as needed
```
### Adding Connections
Define connections between dendrites and the soma:

```matlab
% Add connections using dendrite IDs
% Syntax: model = model.addConnection(ID, List_Proximal, List_Distal);
% Connect dendrite 1 to the soma (ID = 0)
model = model.addConnection(1, 0, []);

% Alternatively, use a connection string
% Example: '0-[1,2]-3' connects the soma to dendrites 1 and 2, which connect to dendrite 3
model = model.addConnectionStr('0-1');
```
### Adding Stimuli
Apply external stimuli to dendrites:

```matlab
% Define a stimulus function (e.g., a sine wave)
StimulusFunction = @(t) sin(2 * pi * 1 * t); % 1 Hz sine wave

% Add the stimulus to dendrite 1
model = model.addStimuli(1, StimulusFunction);
```

### Running Simulations
Prepare and run the simulation:

```matlab
% Set the integration step size
h = 0.01; % Example step size in milliseconds
model = model.buildModel(h);

% Set the maximum simulation time
tmax = 1000; % in milliseconds

% Run the simulation using the RK4 method with the MEX function
[t, solution] = model.RK4CPP(tmax);
```
### Visualization
Visualize the Dendritic Tree
```matlab
% Visualize the dendritic tree in 3D
model.visualize3D();
```
### Plot Simulation Results
```matlab
% Extract membrane potentials (u) from the solution
numDendrites = model.numDendrites;
u = solution(1:2:end, :); % Rows corresponding to 'u' variables

% Plot membrane potentials over time for each dendrite
figure;
for idx = 1:numDendrites
    plot(t, u(idx, :));
    hold on;
end
xlabel('Time (ms)');
ylabel('Membrane Potential (u)');
title('Membrane Potentials Over Time');
legend(arrayfun(@(d) sprintf('Dendrite %d', d.ID), model.dendrites, 'UniformOutput', false));
```
## Project Structure
- ``mexRK4Solver.cpp``: C++ MEX function implementing the RK4 solver.
- ``Dendrite.m``: MATLAB class representing a dendrite.
- ``DendriteTreeModel.m``: MATLAB class for building and simulating the dendritic tree model.
- ``soma_rhs.m``: Function defining the RHS of the ODEs for the soma.
- ``dendrite_rhs.m``: Function defining the RHS of the ODEs for the dendrites.
- ``README.md``: Project documentation (this file).

## Changelog

### [??.??.??] - Next update incoming soon
- **New**:
- **Updated**:
- **Fixed**:
- **Removed**:

---

### [29.09.24] - Finally on git

## Future Plans
- Restructure project:
  ```
  main
  |-src
  | |-matlab
  | | |-*.m
  | |-cpp
  |   |-*.cpp
  |-examples
  | |-*.m examples, covering full functionality
  |-compiled
  | |-*.mexw64
  |-readme
  ```
- Soma-only usage (validate if everything works with current setup, probably some bulky fixes)
  - Optimize
- Results interpreter (probably 'll make some ``.m`` scripts as a starting point
  - Save results in some model-associated way
  - Methods for different plots
- Optional log from cpp
- Snapshots
  - Save model & simulation results in some convenient way (probably after Results interpreter is implemented)
- Check is mex's are compiled on init
- A guide on how to change neuron model
  - Probably a massive overhaul of the project & built-in different models
- Some other cool stuff
