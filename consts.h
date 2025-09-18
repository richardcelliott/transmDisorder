#ifndef _CONSTS_
#define _CONSTS_

#include <iostream>

//Simulation basics.
#define NR 120 //number of dyes in the transmission line, or collocation points.
#define NT 300        //number of discretizations in time. This is
                     //also (naturally) equal to the Ns in Glenn's book.
//Program parameters: frequency of file writing, etc.
#define ENRGY_PTS 1000	//Number of points in time to record the energies.
#define REGPRINT 500 	//Print to file interval for probabilities.
//Energy values for the excitons.
#define ALFAA 1.0     //The ratio wJA/w0, for hopping down line A.
#define ALFAB 1.0     //The ratio wJB/w0, for hopping down line B.
//Shape of the well, energy landscape of excitons.
#define R1 55		//First node in the intrX region.
#define RS 4		//Number of middle interaction sites.
//Relevant initial condition settings. These are for a Gaussian.
#define SIGMA 6.0     	//Gaussian width in lattice spacings.
#define XX 5     	//Width of Bound region in width sigma
#define KR 30           //Momentum 'mode' for initial conditions in the R-line.



extern double pie;

#endif //  _EXCITON_PROP_
