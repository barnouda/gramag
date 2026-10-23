# GraMag - Gravimetry and Magnetism modelling

## Overview

These codes compute the vertical gravity field or the magnetic field along a given direction for a layer model surrounding a grid model.

The four versions of the code:
* gralite: Computes the vertical gravity field for a layer model surrounding a grid model.
* grafull: Does the same and additionnally provides the linear coefficients that allows to recompute the gravity field with different values of density at the grid nodes.
* maglite : Computes the magnetic field in a given direction for a model of magnetizations in a given direction.
* magfull : Does the same and additionnally provides the linear coefficients that allows to recompute the magnetic field with different values of magnetization at the grid nodes.

## Authors

Created by Olivier Coutant in 2010-2011.

Modified by Olivier Coutant and Anne Barnoud in 2013-2015.

## Where is the code?

The GraMag code is available on Codeberg: https://codeberg.org/barnouda/gramag

## References

When using this code, please cite:
* Coutant O., Bernard M.-L., Beauducel F., Nicollin F., Bouin M.-P., Roussel S., Joint inversion of P-wave velocity and density, application to La Soufrière of Guadeloupe hydrothermal system, Geophysical Journal International, Volume 191, Issue 2, November 2012, Pages 723–742, https://doi.org/10.1111/j.1365-246X.2012.05644.x

and one of the following two citations depending on your use:
* [for gravimetry] Barnoud A., Coutant O., Bouligand C., Gunawan H., Deroussi S., 3-D linear inversion of gravity data: method and application to Basse-Terre volcanic island, Guadeloupe, Lesser Antilles, Geophysical Journal International, Volume 205, Issue 1, April 2016, Pages 562–574, https://doi.org/10.1093/gji/ggw030
* [for magnetism] Barnoud A., Bouligand C., Coutant O., Carlut J., Magnetic structure of Basse-Terre volcanic island (Guadeloupe, Lesser Antilles) inferred from 3D inversion of aeromagnetic data, Journal of Volcanology and Geothermal Research, Volume 348, December 2017, Pages 1-11, https://doi.org/10.1016/j.jvolgeores.2017.10.005
