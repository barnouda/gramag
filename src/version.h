/*
 *  Versions prog GraMag
 *
 *  Created by Olivier Coutant on 04/03/11.
 *  Université Joseph-fourier, Grenoble.
 *
 *
 *  V1.1  09/04/13  OC	Memory savings : X, Y and DXY float and uchar instead of double.
 *
 *  V1.2  09/04/13  AB	Debug V1.1 : values of discretisation steps in DXY corrected
 *  			New option -nc : if present, do not check if receiver is located inside the topography boundaries
 *  	  17/04/13  AB	Same modifications for grafull
 *
 *  V1.3  21/05/13  OC  Lines dx, dy and lat, lon removed from parameters file
 *                      The input grid is now in volumic mass (instead of velocity)
 *                      Rq: so rhoVp.c is not used anymore (previously used to convert velocity model to volumic mass)
 *
 *  V1.4  21/05/13  AB  Bug corrected in topprismLowRes de topprism9.c (line 460 *rho removed)
 *
 *  V1.5  26/08/13  OC  Modify call of input file
 *                      Bugs corrected in the case where no topography is given
 *                      Bug corrected in primJ1 when a=0 (call atan2 rather than atan)
 *                      Bug corrected in the case where the grid is completly beneath the flat topography (prismIntegration)
 *
 *  V1.6  21/01/14  AB  Magnetism: maglite et magfull (adaptation of the code to magnetism field computation)
 *
 *  V1.7  12/02/15  OC  New test on the size of the grid
 *
 *  V1.8  13/08/15  OC  Bug corrected : add a test in loop8/9 to shift the receiver when it is located at a corner or on the side of a prism
 *
 *  V1.9  01/08/25  AB  Ten years later, clean code for github
 *
 */
#define VERSION "1.9 01/08/2025"
