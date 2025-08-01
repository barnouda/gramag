# GraMag - Gravimetry and Magnetism modelling

These codes compute the vertical gravity field or the magnetic field along a given direction for a layer model surrounding a grid model.

## Authors

Created by Olivier Coutant in 2010/2011.

Modified by Olivier Coutant and Anne Barnoud in 2013-2015.

## References

When using this code, please cite:
* Coutant O., Bernard M.-L., Beauducel F., Nicollin F., Bouin M.-P., Roussel S., Joint inversion of P-wave velocity and density, application to La Soufrière of Guadeloupe hydrothermal system, Geophysical Journal International, Volume 191, Issue 2, November 2012, Pages 723–742, https://doi.org/10.1111/j.1365-246X.2012.05644.x

and one of the following two citations depending on your use:
* [for gravimetry] Barnoud A., Coutant O., Bouligand C., Gunawan H., Deroussi S., 3-D linear inversion of gravity data: method and application to Basse-Terre volcanic island, Guadeloupe, Lesser Antilles, Geophysical Journal International, Volume 205, Issue 1, April 2016, Pages 562–574, https://doi.org/10.1093/gji/ggw030
* [for magnetism] Barnoud A., Bouligand C., Coutant O., Carlut J., Magnetic structure of Basse-Terre volcanic island (Guadeloupe, Lesser Antilles) inferred from 3D inversion of aeromagnetic data, Journal of Volcanology and Geothermal Research, Volume 348, December 2017, Pages 1-11, https://doi.org/10.1016/j.jvolgeores.2017.10.005

## Description of the density/magnetization model

Grid of nodes surrounded by a layered model

```
______________________________________________________

               .  .  .  .  .  .  .
               .  .  .  .  .  .  .     first layer
               .  .  . _._ .  .  .
               .  .  ./ . \.  .  .
___topography__.__.__/  .  \__.__.__        __________
               .  .  .  .  .  .  .  \______/
               .  .  .  .  .  .  .  
_______________.  .  .  .  .  .  .____________________
               .  .  .  .  .  .  .
               .  .  .  .  .  .  .     second layer
               .  .  .  .  .  .  .
_______________.  .  .  .  .  .  .____________________

                                       third layer
                                       (half-space)

```

Layers have constant densities or magnetizations outside the grid.

The grid model can be initialized from a density/magnetization model or it can be initialized to the same values as the extern layered model.

Within the grid, densities/magnetizations are defined everywhere as the trilinear interpolation of values assigned at the nodes.

## Repository content

* code sources in src
* makefile to build src in build and bin
* example with a toy cone synthetic example

## Codes

* gralite: Computes the vertical gravity field for a layer model surrounding a grid model.

* grafull: Does the same and additionnally provides the linear coefficients that allows to recompute the gravity field with different values of density at the grid nodes.

* maglite : Computes the magnetic field in a given direction for a model of magnetizations.

* magfull : Does the same and additionnally provides the linear coefficients that allows to recompute the magnetic field with different values of magnetization at the grid nodes.

## Common options

* -l: If present, assign the outside layer density values to the grid nodes.

* -np: If present, do not project the receiver altitude on the topography. This allows to compute the gravity/magnetic field for receiver located above the topography (e.g. for airborne data).

* -nc: If present, do not check if receiver is located inside the topography boundaries. This allows to compute the gravity field of a partial topography model away from the receiver (e.g. the effect of bathymetry away from the network).


## Input Files

Here is a short description of the input files. Examples in the example folder can be useful as well.

* The input parameter file contains the following lines:
	1. prefix (prefix name used for output files)
	2. topo.bin (topography file using the topo format given below, or "none" if no topography given)
	3. model.par (grid model parameters)
	4. model.bin (velocity grid model using the format given below, or "none" if no velocity given)
	5. data.txt (receiver locations where we compute the vertical gravity field. (x,y,z) in ascii, z is altitude, positive upward)
	6. 630922 654922 (x integration limit for density model, in meter, (e.g, UTM coordinates))
	7. 1763803 1794803 (y integration limit for density model, in meter, (e.g, UTM coordinates))
	8. -2000 0 (zmin and zmax depth of integration; positive downward, origin at geoid. The upper limit is either given by zmin or by the topography, the lower limit zmax must of course be located below the lowest topography point)
	9. 2 (number of layers for the layered density model outside the grid)
	10. 2160 -1050 (lines 10 and following, for every layers: volumic mass and bottom depth)
	11. 2710 0

* Topography binary file:
	* binary file written in C or fortran option open (.... FORM = 'unformatted', ACCESS = 'stream')
	* In our case, x is longitude in UTM coordinates in meters
	
	containing:

	* ymin  (float or real*4) 
	* ymax  (float or real*4)
	* xmin	(float or real*4)
	* xmax	(float or real*4)
	* step_x_and_y	(float or real*4)
	* nx	  (int or integer*4)
	* ny	  (int or integer*4)
	* topo	(nx\*ny float or real\*4) (fast direction is x, i.e. longitude first if x is long)

* Model parameter text file with the following lines:
	1. 1 (because it's like this)
	2. 642600.000000 1773800.000000 -1500.000000 (xmin ymin zmin, z downwards)
	3. 15  15  11 (number of nodes in the x, y, z directions)
	4. 100.000000  100.000000  100.000000 (xstep ystep zstep)

* Model binary file with densities/magnetizations at nodes (float)

* Text data file containing the x, y, z locations of the data (z upwards!)

* For magnetism (maglite/magfull), "dir.par" text file (fixed name) containing two lines:
	1. xnorthward yeastward zdownward of rock magnetization
	2. xnorthward yeastward zdownward of the regional magnetic field at the time of observations

## Usage

* Prepare required input files (see example)
* Adapt makefile if needed
* Build code (e.g. make gralite)
* Run code (e.g. ./bin/gralite -l -np -nc input_parameters.in)
* If using these codes, please cite the above references.

## Example (./example)

The example dataset is a synthetic volcano (conic topography) with laterally varying densities/magnetizations.

The synthetic topography, model, data and input files for modelling are generated with the example_toy_cone.m matlab code and its dependencies. You may ignore the *muog* files for use with another code.

To run the code with the example dataset :
* build code with the makefile
* run the code in the example folder, ex:
	* ../bin/grafull -np -nc input_gravi_modelling.in
	* ../bin/maglite -np -nc input_mag_modelling.in

