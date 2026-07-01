## Phase field fracture, fatigue, and fatigue simulation acceleration models

This repository implements a phase field fracture & fatigue model based on the finite element library [deal.II](https://www.dealii.org/), and six fatigue simulation acceleration models introduced in our paper (UPDATE AFTER SUBMISSION). 

The repository features:

* Easy implementation of multiphysics problems: For a new physics added to the existing system, one only needs to implement (1) how residual (or the Right Hand Side) and tangent stiffness matrix (or stiffness matrix) are calculated on quadrature points, (2) calculating visualized fields, and (3) defining the staggered scheme.

* Easy problem set-up: Problem definition with (1) an ABAQUS .inp file defining the mesh and boundaries, (2) a text file describing Dirichlet or Neumann boundary conditions (based on ``Surface'' sets defined in Abaqus .inp file), and (3) a parameter file for computing setups, preferably without modifying any code.

* Adaptive mesh based on gradients of the phase field. 

* Parallelism with MPI, working smoothly with multiple nodes on high-performance computing clusters. 

* Extensible interfaces for new coupling field and fatigue simulation acceleration algorithms. 

### Installation

We recommand using `docker` or `singularity` containers to avoid dependency issues when installing `deal.II`. This code is tested on the official deal.II docker image [v9.6.0-noble-amd](https://hub.docker.com/layers/dealii/dealii/v9.6.0-noble-amd64/images/sha256-b69efa3a16499913927e062b71104d2e8e6780a7fe387d91008b87c1eda64334) and [v9.4.0-focal](https://hub.docker.com/layers/dealii/dealii/v9.4.0-focal/images/sha256-ca02d4f9e2cc6fbae5226f3c6fb4c6ff9399770ec9c376f2bbec893139341367).

### Quick start

This command compiles the code and runs a single notch tension test (comparable to Miehe et al., CMAME, 2010) with 8 MPI threads:

```shell
./compile_and_run.sh -n 8 -f parameters/singleNotchTension.prm
```

### Reproducing the results in the paper

`UPDATE_AFTER_SUBMISSION.sh` contains all the commands to reproduce the results in the paper. It is recommended to run the commands one by one, as some of them are time-consuming.

### Usage

* Mesh file: It uses meshes defined in ABAQUS-generated .inp files (see `meshes/singleNotchDense.inp` for example).
* Boundary conditions: First, define "Surfaces" in ABAQUS (in the Assembly module). Then define a boundary configuration file (see `meshes/singleNotchTension_boundary.txt` for example). The configuration is defined according to a specific format:

	* Each line defines: Surface-ID, type of constraint, constrained dof, and value(s).
	* For Dirichlet boundaries:
	  * For `velocity`, the fourth number is in mm/s (or the unit of the rate of the field variable)
	  * For `dirichlet`, the fourth number is in mm (or the unit of the field variable)
	  * For `sinedirichlet`/`triangulardirichlet`, the fourth part is frequency(Hz), mean(mm), and amplitude(mm), for example "20 1 2"
	* For Neumann boundaries: the third part is a series of floats denoting the vector of the neumann boundary (set 0 if is a scalar field).
	  * For `neumann`, the third part is in MPa (or the unit of gradient)
	  * For `neumannrate`, the third part is in MPa/s (or the unit of the rate of gradient)
	  * For `sineneumann`/`triangularneumann`, the third part is dimensionless. And add a fourth part being frequency (Hz), mean (MPa), and amplitude (MPa), for example "20 1 2"
* Parameters: See `parameters/singleNotchTension.prm` for an example of parameter definitions. All available parameters are shown in `include/parameters.h`. These parameters are straightforwardly named.

* Compile the code and execute a project using a parameter file:

  ```shell
  ./compile_and_run.sh -n [No_MPI_THREADS] -f [PARAMETER_FILE]
  ```

  Some additional arguments are available:

  *  `-s`:  `true` or  `false`, meaning whether to use the previously compiled code.
  *  `-r`:  `debug` or  `release`. 

### Extensions

If you want to define a new field, ideally you have to (and only have to) imitate  `phase_field.h` or  `elasticity.h` (that inherits  `abstract_field.h`) to define operations for the field, and then imitate  `phase_field_fracture.h` (that inherits `abstract_multiphysics.h`) to define a regime of solving the multiphysical system (particularly the staggered scheme). 

If you want to define a new fatigue simulation acceleration algorithm, see `include/adaptive_timestep.h` and `include/fatigue_degradation.h` for our implementations as references.

### Cite us

If you use this code in your research, please cite our paper (UPDATE AFTER SUBMISSION).
