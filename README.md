## deal.II phase field fracture

### Install

This code is tested on the official deal.II [docker image (v9.4.0-focal)](https://hub.docker.com/r/dealii/dealii/tags).

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
  ./compile_and_run.sh -n 8 -f parameters/singleNotchTension.prm
  ```

  where  `-n` defines the number of MPI processes and  `-f` defines the parameter file

  There are some additional parameters:

  *  `-s`:  `true` or  `false`, meaning whether to use the previously compiled code.
  *  `-r`:  `debug` or  `release`. 

### Extensions

If you want to define a new field, ideally you have to (and only have to) imitate  `phase_field.h` or  `elasticity.h` (that inherits  `abstract_field.h`) to define operations for the field, and then imitate  `phase_field_fracture.h` (that inherits `abstract_multiphysics.h`) to define a regime of solving the multiphysical system (particularly the staggered scheme). 

 
