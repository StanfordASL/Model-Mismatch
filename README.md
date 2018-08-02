# Model-Mismatch

Code to accompany: Robust Tracking with Model Mismatch for Fast and Safe Planning: an SOS Optimization Approach

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. All code is written in MATLAB.

### Prerequisites

There are three packages required: 
* [spotless](https://github.com/spot-toolbox/spotless) : polynomial optimization parser, 
* [MPT 3.0](https://www.mpt3.org/) : manipulating and plotting convex shapes,
* [Mosek](https://www.mosek.com/) : SDP solver for the SOS programs.

All packages are freely available at the provided links.

### Installing

Having installed the prerequisites (and adding them to the MATLAB path), download the repo, navigate to the repo in MATLAB, and execute the following in the command window:

```
run startup_mm.m
```

This will add all necessary sub-directories of the repo to the MATLAB path.

## Code Structure

The code is divided into a core computational backend and custom code for each example.

Backend:
* optimizers/ : contains the SOS problem definitions for the K and V subproblems (using the spotless toolbox),
* auxiliary/ : contains auxiliary functions such as converting msspoly (spotless) functions to MATLAB anonymized functions, plotting ellipsoids, computing Chebyshev polynomial approximations for non-polynomial functions.

The custom code for each example (Dubins, Car5D, PVTOL, Plane8D) contains the problem specific definitions such as dynamics, state and control constraints, polynomial degrees, etc. As an example, the following highlights the key components of each of these example folders. This should serve as a template to setup your own tracking/planning system. 

Consider the 5D car vs 3D Dubins example in the paper (section 5.1), corresponding to the Car5D folder. The main run file is *car5D.m*. The structure of this file is:

* tracking control constraints: lines 3-17
* setup indeterminates for SOS program: lines 19-25
* dynamics definition, including any Chebyshev function approximations: lines 27-49
* bounded state and Jacobian: lines 51-64
* state and planning control constraints: lines 66-77
* setup tolerances and degrees: lines 79-90
* run initialization algorithm: line 95
* main algorithm loop (Algorithms 2 and 3): lines 114-178
* numerical check on solution: line 183
* save solution for simulation: lines 189-199
* plots and projections, line 200 onwards. 

The initialization function (*initialize_car5D.m*) runs the infeasible start algorithm (Algorithm 1) and calls *solve_car5D_start_1.m* and *solve_car5D_start_2.m* to run the K and V sub-problems. The main algorithm calls *solve_car5D_1.m* and *solve_car5D_2.m* to run the K and V sub-problems. These functions provide the interface between the SOS programs, the solver MOSEK, and MATLAB. 

### Simulation and Verification

The simulation folder in each example folder is labelled with "_Sim." The main run file will be labelled as *test_<system>.m*. The file loads the V and K functions, generates or loads a planned path for the planning system, simulates the system using the tracking dynamics and the controller K, and generates evaluation plots. 

Datafiles to replicate all results in the paper are provided within the repo, except for the FMT* graph structure which will be (deterministically) generated upon first execution of *test_Plane8D.m* in the Plane8D/Plane_Sim folder. 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

