#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <complex>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <armadillo>
#include "consts.h"
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
using namespace arma;


#define qudi(x) #x
#define numberize(x) qudi(x)

vec timeEvolv(string, string, string, cx_vec, vec, string, double, uword);
cx_vec create_init_conds();
cx_vec create_init_conds_from_input(string,double,double);
cx_vec create_static_init_conds(double,double);
bool hasEnding(std::string const &, std::string const &); 
vec create_diag_elements(double);
double get_norm_of_1d();
void probs_to_file(string,int, cx_vec);
void flucts_to_file(string, vec, uword);
void wavefunc_to_file(string,int,cx_vec);
double correl8(string,int,cx_vec,cx_vec);
void dump_aves2file(string,cx_vec,string);
void dump_aves2file_real(string,vec,string);
cx_vec create_bound_initConds(double);
cx_vec create_bound_initConds2(double);
std::string create_output_dir(double,double,string);
std::string create_subdir(string, string, uword);
std::string create_subdir_gen(string, string);

double initTime = 0.00;
double pie = 2.0 * acos(0.0);
complex<double> ii = complex<double>(0,1.0);
uword r2 = R1+RS;
double r0 = 109.5;	//make a half integer for an FFT used below (NR even).
double endTime = 42.0;
double ds = (1.*endTime-initTime)/(1.0*(NT-1)); 

int main(int argc, char* argv[])
{

    //number of statistical samples to use for each noise window.
    uword nW = 1111;

    //Type of simulation: noise in site energy (E) or in hopping (J).
    //string typ_ = "J";
    string typ_ = "E";

    uword nMag = 32;
    double stepSze = 0.04;
    //loop over the noise magnitude size.
    for(uword kk=0; kk<nMag; ++kk) {

	//Size of the window for noise in the middle section.
	double initDelta_ = 0.01+1.0*kk*stepSze;

	//create the main results directory.
	string mainDir = create_output_dir(initDelta_,1.*KR,typ_);
	//create another secondary results subdirectory.
	string subDir01 = "flucts";
	string fluctsSubDir = create_subdir_gen(mainDir, subDir01);
	//create another secondary results subdirectory.

	string subDir03 = "probs";
	string probsSubDir= create_subdir_gen(mainDir, subDir03);

	vec arrFluct(NR, fill::zeros); 
	vec corrVec(NT, fill::zeros); 
	vec corrVecAve(NT, fill::zeros); 
	//mat fluctsMatrx(RS, nW, fill:zeros);

	//Main statistical loop.
	for(uword i=0; i<nW; ++i) {
	    arma_rng::set_seed_random();
	    //vec flucts( RS, fill::randn );	//variance is \sigma^2=1 here.
	    vec flucts = randn(RS, distr_param(0.00,initDelta_)); 
	    			//Here, stdev=initDelta, ave=0.
	    //fluctsMatrx.col(i) = 1.0 + initDelta_*flucts;
	    for(uword k=0; k<R1; ++k) { arrFluct(k) = 3.0; }
	    for(uword k=R1; k<r2; ++k) { 
	        arrFluct(k) = 3.0 + flucts(k-R1); 
	        //arrFluct(k) = 1.0;
	        //cout << arrFluct(k) << endl;
	    }
	    for(uword k=r2; k<NR; ++k) { arrFluct(k) = 3.0; }
	    flucts_to_file(fluctsSubDir, arrFluct, i);

	    //create a secondary results subdirectory. One for each iter/seed.
	    //
	    string subDir02 = "wavefunc";
	    string wavefuncSubDir = create_subdir(mainDir, subDir02, i);
	    //cout << probsSubDir << endl;

            //initial conditions.
	    cx_vec psi0 = create_init_conds();
	    //Dump it to file.
	    //probs_to_file(probsSubDir,0,psi0);
	
	    //create another secondary results subdirectory.
	    //string subDir02 = "correl";
	    string correlSubDir = "junk";
	    //string correlSubDir = create_subdir_gen(mainDir, subDir02);
	    //cout << correlSubDir << endl;

	    //double grbge = correl8(correlSubDir,0,psi0,psi0);
	    //cout << endl;
	    //cout << grbge << endl;
	    //cout << endl;

	    //This returns the maximum of the cross-correlation vector as a
	    //function of time, for each statistical sample.
	    corrVec = timeEvolv(probsSubDir,correlSubDir,wavefuncSubDir,psi0,arrFluct,typ_,initDelta_,i);
	    //Sum them: an average will be recorded in another file.
	    corrVecAve += corrVec;
	}//closes statistical loop.

	//Divide by the size of the statistical sample.
	corrVecAve /= (1.0*nW);

	//Dump average to file. Not necessary, nor needed, but here it is.
	//string aveDum = "correl8t_ave";
	//dump_aves2file_real(fluctsSubDir,corrVecAve,aveDum);

    }//closes the loop over noise window size.

	//corrVecAve.print("corrAve");


	return 0;

}

vec timeEvolv(string probsDir, string correlDir, string wavefuncDir, cx_vec psiInit, vec arrFluct, string typ, double initDelta, uword i)
{

	//Coefficients in the matrix solver.
	//double b = (1.*ALFAA);
	//double c = (1.*ALFAB);

	/*
	cout <<  endl;
	cout << "    ===================================================" << endl;
	cout <<  endl;
	cout << "        Run has begun. Type:  "<< typ << "\t\t" << "Parameters: " << endl;
	cout <<  endl;
	cout << "            En. Delta: " << endl;
        cout << "                          " <<	initDelta << endl; 
	cout << "            HopRight: " << endl;
        cout << "                          " <<	(1.*ALFAA) << endl; 
	cout << "            HopLeft: " << endl;
        cout << "                          " <<	(1.*ALFAB) << endl; 
	cout << "            NR-1, R1, R2: " << endl;
	cout << "                          " << NR-1 << ", " << R1 << ", " << r2 << endl; 
	cout << "            R0, KR-1: " << endl;
	cout << "                          <==" << r0 << ", (" << KR-1 << "), " << endl; 
	cout <<  endl;
	cout << "    ===================================================" << endl;
	cout <<  endl;
	*/

        //Create the Ham matrix first. This is all STATIC. Note.
	mat tDiagMatrx(NR, NR, fill::zeros);


	if( typ == "J" ) {
		//cout << typ << endl;

		vec arrFluctMinusOne(NR-1,fill::zeros);
		//all values except the last, which ends up on the corners.
		for(uword k=0; k<arrFluct.n_elem-1; ++k) {
			arrFluctMinusOne(k)  = arrFluct(k);
		}
		//double aFluctEnd = arrFluct.tail(1);

		tDiagMatrx.diag(1) = arrFluctMinusOne;	//offdiagonals.
		tDiagMatrx.diag(-1) = arrFluctMinusOne;	//offdiagonals.
		tDiagMatrx.diag(-(NR-1)) += arrFluct.tail(1); //corner values.
		tDiagMatrx.diag(NR-1) += arrFluct.tail(1);	//corner values.
		//tDiagMatrx.print("tDiag:");

		//Set diagonal here.
		vec enDiag(NR, fill::zeros); 
		for(uword k=0; k<NR; ++k) { enDiag(k) = 1.0; }
		//Create the vector diagonal. Somewhat complicated here.
        	//vec vDiag = create_diag_elements();
		tDiagMatrx.diag() = enDiag;


	}
	else if( typ == "E" ) {
		//cout << typ << endl;

		vec arrFluctMinusOne(NR-1,fill::zeros);
		//all values except the last, which ends up on the corners.
		for(uword k=0; k<arrFluct.n_elem-1; ++k) {
			arrFluctMinusOne(k)  = arrFluct(k);
		}

		tDiagMatrx.diag(1) += 1.0;	//offdiagonals.
		tDiagMatrx.diag(-1) += 1.0;	//offdiagonals.
		tDiagMatrx.diag(-(NR-1)) += 1.0; //corner values.
		tDiagMatrx.diag(NR-1) += 1.0;	//corner values.
		//tDiagMatrx.print("tDiag:");

		//Set diagonal here.
		vec enDiag(NR, fill::zeros); 
		for(uword k=0; k<NR; ++k) { enDiag(k) = arrFluct(k); }
		//Create the vector diagonal. Somewhat complicated here.
        	//vec vDiag = create_diag_elements();
		tDiagMatrx.diag() = enDiag;



	}
	else {
		cout << "bork."<< endl;
	}



	//Timer.
	wall_clock timer;

	cout << endl;
	cout << "====================================" << endl;
	cout << "  Starting matrix diagonalization."  << endl;
	cout << "====================================" << endl;
	timer.tic();

	vec eigval;
	mat eigvec;

	eig_sym(eigval, eigvec, tDiagMatrx);

	double nSec = timer.toc();
	cout << endl;
	cout << "========================================" << endl;
	cout << "Matrix diagonalized in " << nSec << " seconds. This is " 
			<< nSec/60.0 << " minutes." << endl;
        cout <<	"    Now starting time evolution..."  << endl;
	cout << "========================================" << endl;

	cout << endl;
	cout << endl;

	cx_vec cN(NR, fill::zeros); 
	//cout << "size of eigvec: " << eigvec.n_cols << endl;
	//cout << "col of eigvec: " << endl;
	//eigvec.col(12).brief_print("12th vec stuff");
	//exit(1);

	//Take initial conditions and decompose into elements,
	//This determines the coefficients cN.
	for(uword k=0; k<eigvec.n_cols; ++k) {
		vec z = eigvec.col(k); // extract a column vector
    		cN(k) = dot(psiInit,z);
		//cout << k << "\t" << z << endl;
	}
	//cN.print("Coeff");

	vec corr8_t(NT, fill::zeros); 
	cx_vec energyHop(NT, fill::zeros); 
	cx_vec energyActv8A(NT, fill::zeros); 
	cx_vec energyTot(NT, fill::zeros);

	cx_vec psiRtSav(NR, fill::zeros); 

	wall_clock timer2;
	cx_vec psiUpd(NR, fill::zeros); 
	double tMe;
	timer2.tic();

	//Next, the time iteration.
	cout << "================================================" << endl;
	cout << "Next, eigenvector manipulations and file dump." << endl;
        for(int j = 1; j < NT; j++)
        {
		tMe = ds*j;
		for(uword k=0; k<eigvec.n_cols; ++k) {
			double eN = eigval(k);
			//cout << "eigs: " << eN << endl;
			psiUpd(k) = cN(k)*exp(-1.0*ii*eN*tMe); 
					//each k is scalar.
		}

		//Here we create the vector of wavefunction values, psiRSt,
		//for the full wavefunction.
		cx_vec psiRt(NR, fill::zeros); 
		for(uword k=0; k<eigvec.n_cols; ++k) {
			vec z = eigvec.col(k); // extract a column vector
			psiRt += psiUpd(k)*z;	
		}



		//Dump probs to file for whole wavefunc.
		//probs_to_file(probsDir,j,psiRt);

		//Correlation calculation.
		//correl8(correlDir,j,psiInit,psiRt);


		//corr8_t(j) = correl8(correlDir,j,psiInit,psiRt);
		psiRtSav = psiRt;
	} //closes loop over time evolution (indexed by j).

	wavefunc_to_file(wavefuncDir, NT-1, psiRtSav);


	//string sT = to_string(i);
	//string aveDum = "correl8t_" + sT;
	//dump_aves2file_real(correlDir,corr8_t,aveDum);




	double nSec2 = timer2.toc();
	cout << "    File dump finished in  " << nSec2 << " seconds. This is " 
			<< nSec2/60.0 << " minutes." << endl;
	cout << "================================================" << endl;
	cout << endl;



	cout << endl;
	return corr8_t;
}

vec create_diag_elements(double w21val)
{

	vec wR(NR, fill::zeros);
	vec wS(NR, fill::zeros);
	mat ham(NR*NR, NR*NR, fill::zeros);

	//cout << R1 << "\t" << r2 << endl;
	//exit(1);
	//The limited range with R1 < r < R1+2 having different values.
	//Outside the well, to the left.
	for(uword k=0; k<R1; ++k) {
		wR(k) += 0.;
		wS(k) += 0.;
	}
	//Inside the well. Last node is R2. Not less than R2.
	for(uword k=R1; k<r2; ++k) {
		wR(k) += (1.0*w21val);
		wS(k) += (1.0*w21val);
	}
	//Outside the well, to the right. Starts at R2+1, goes to NR-1.
	for(uword k=r2; k<NR; ++k) {
		wR(k) += 0.;
		wS(k) += 0.;
	}
	

	for(uword i=0; i<(NR*NR); ++i) {
		//indices from the grid assignment: i = r * N_r + s.
		uword s = i % NR;	//'fast' index.
		uword r = (i-s)/NR;	//slow, stride index.

		for(uword j=0; j<(NR*NR); ++j) {
		uword sprime = j % NR;;
		uword rprime = (j-sprime)/NR;
		
		if( (s == sprime) && (r==rprime)) {
			ham(i,j) += ( wR(r)+wS(s) );

	//ADD INTERACTION TERMS. These are weird, indexed at some diagonals.
			//if( r==s) {
			//	ham(i,j) += (-1.*BTA);
			//}
		}

		}
	}
        vec diags = diagvec(ham);
	//diags.print("diag");
	//
	//for(uword i=0; i<(NR*NR); ++i) {
		//cout << i << "\t" << diags(i) << endl;
	//}
	//exit(1);

	return diags;
}

cx_vec create_init_conds()
{

	//Declare a matrix, fill with IC profile, then convert to a vector.
	cx_vec icV(NR, fill::zeros);

	//Create a Gaussian at (r0,s0), and symmetrize the grid about this
	//location. We must make this periodic in the grid. this is not easy.
	
	//First create a Gaussian of the proper widths in (r,s) at the MIDPOINT
	//of the grid, then translate to the proper peak location w/FFT.  
	double midX = (1.*(NR-1)/2.0);
        for(uword r = 0; r < NR; r++) { 
	  icV(r) = exp(-1.0*((1.*r-midX)*(1.*r-midX))/(2.*SIGMA*SIGMA));
	}

	//Normalize the square of a 1-dim wavefunc, the symmetric Gaussian. A
	//little dirty here.
	double nrmSq = get_norm_of_1d();

	//Fourier Transform.
	cx_vec icV_k = fft(icV);

	//Multiply by a plane wave for translation of the Gaussian peak to
	//(r0). Translation values first:
	//
	double nR_trnsl = 1.0*midX-1.0*r0; 
	//cout << nR_trnsl << endl;
	//(1.*midX+5.0); //1.*(1.*midX-1.*r0);
	//Multiply by a plane wave with the translation vector.
	//Phase sign is correct by trial and error (not sure of FFT convention
	//with armadillo, it is not well-documented).
	for(uword kr=0;kr<NR;kr++)
	{
	  double kr_ = 1.0*kr;
	  icV_k(kr) *= exp(1.0*ii*((2.0*pie*kr_*nR_trnsl)/(1.0*NR)));
	}

	//Now transform back to real space for a Gaussian at the proper
	//location of the initial condition.
	cx_vec icM = ifft(icV_k);

	//Multiply by the plane wave in both dims. And also normalize.
	for(uword r=0;r<NR;r++)
	{
	  icM(r) *= exp(1.0*ii*((2.0*pie*KR*r)/(1.0*NR)));
	  icM(r) /= sqrt(nrmSq);
	}
	//icM.print("icM_Init:");

	return icM;

}

bool hasEnding (std::string const &fullString, std::string const &ending) 
{

    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } 
    else {
        return false;
    }

}

double correl8(string sDir, int iT, cx_vec psi0, cx_vec psiT)
{

        //const char *fileNameC;
        //ofstream fileStreamC;

	double maxCorr;
	//double tiMe = initTime + iT*ds;
	//string sT = to_string(tiMe);

	//std::string s0c = sDir + "/" + "correl_" + sT + ".dat";
	//fileNameC = s0c.c_str();

        //fileStreamC.open(fileNameC);
        //if ( !fileStreamC.is_open() )
        //{
                //cout << "fileStream for C is broken" << endl;
                //cout << endl;
                //exit(1);
        //}
        //else
        //{


	    cx_vec corr(NR, fill::zeros);
	    cx_vec nrm1(1, fill::zeros);

            for(uword v = 0; v < NR; ++v) { 
                    for(uword s = 0; s < NR; ++s) { 
	                uword t=s+v;
			//reflect through BCs.
			t = t - floor(t / NR) * NR;
			//cout << v << "\t" << s << "\t" << t << endl;
	                corr(v) += conj(psi0(s))*psiT(t);    
	                //corr(v) += (psi0(s)*conj(psi0(s)))*(psiT(t)*conj(psiT(t)));    
		    }
	            nrm1 += conj(psi0(v))*psi0(v);    
	    }
	//    cout << corr(0) << "\t" << nrm1 << endl;
	 //   exit(1);

	    vec abs_corr(NR, fill::zeros);
	    abs_corr = abs(corr);

	    double nrm = 0.0;
            for(uword j = 0; j < NR; ++j) { 
		    nrm += sqrt( real(corr(j))*real(corr(j)) + imag(corr(j))*imag(corr(j)) );
	    }
	    //cout << nrm << endl;
	    //exit(1);

	    /*
	    //print each contracted probability to file.
	    for(uword j=0; j<corr.n_elem; ++j)
	    {
                //fileStreamC << j << "\t" << real(corr(j)) << endl;
                //fileStreamC << j << "\t" << sqrt( real(corr(j))*real(corr(j)) + imag(corr(j))*imag(corr(j)) ) << endl;
                fileStreamC << j << "\t" << abs_corr(j) << endl;
	    }
	    */
	    //cout << abs_corr.max() << endl;
	    maxCorr = abs_corr.max();

            //fileStreamC << endl;
            //fileStreamC.close();
	//}


	return maxCorr;
}

void dump_aves2file_real(string dirNme, vec dens, string aveDum)
{

        /* char fileNameD[256];
        ofstream fileStreamD;

        sprintf(fileNameD,"%s.dat",aveDum.c_str());

	cout << aveDum.c_str() << endl;
	    for(uword j=1; j<5; ++j) {
		    cout << real(dens(j)) << "\t" << imag(dens(j)) << endl;
	    }
        */

        const char *fileNameD;
        ofstream fileStreamD;

	//double tiMe = initTime + iT*ds;
	//string sT = to_string(tiMe);

	std::string s0c = dirNme + "/" + aveDum + ".dat";
	fileNameD = s0c.c_str();

        fileStreamD.open(fileNameD);
        if ( !fileStreamD.is_open() )
        {
                cout << "fileStream for D is broken" << endl;
                cout << endl;
                exit(1);
        }
        else
        {

	    //print wavefunction correlation at t=0, which is 1 (checked above)
	    if(aveDum == "correl8t_") {
            	fileStreamD << 0.00 << "\t" << 1.0 << endl;
	    }
	    //print each contracted probability to file.
	    for(uword j=1; j<dens.n_elem; ++j)
	    {

	        double tiMe = initTime + j*ds;
                fileStreamD << tiMe << "\t" << dens(j) << endl;
	    }

            fileStreamD << endl;
            fileStreamD.close();
	}

	return;
}

void dump_aves2file(string dirNme, cx_vec dens, string aveDum)
{

        /* char fileNameD[256];
        ofstream fileStreamD;

        sprintf(fileNameD,"%s.dat",aveDum.c_str());

	cout << aveDum.c_str() << endl;
	    for(uword j=1; j<5; ++j) {
		    cout << real(dens(j)) << "\t" << imag(dens(j)) << endl;
	    }
        */

        const char *fileNameD;
        ofstream fileStreamD;

	//double tiMe = initTime + iT*ds;
	//string sT = to_string(tiMe);

	std::string s0c = dirNme + "/" + aveDum + ".dat";
	fileNameD = s0c.c_str();

        fileStreamD.open(fileNameD);
        if ( !fileStreamD.is_open() )
        {
                cout << "fileStream for D is broken" << endl;
                cout << endl;
                exit(1);
        }
        else
        {

	    //print each contracted probability to file.
	    for(uword j=1; j<dens.n_elem; ++j)
	    {

	        double tiMe = initTime + j*ds;
                fileStreamD << tiMe << "\t" << real(dens(j)) << endl;
	    }

            fileStreamD << endl;
            fileStreamD.close();
	}

	return;
}
double get_norm_of_1d()
{

	vec ic1dim(NR, fill::zeros);
	uword midX = NR/2;
	
	//NOTE SQUARE OF GAUSSIAN. OMG.
        for(uword r = 0; r < NR; r++) { 
	ic1dim(r) = exp(-1.0*((r-midX)*(r-midX))/(2.*SIGMA*SIGMA));
	ic1dim(r) *= exp(-1.0*((r-midX)*(r-midX))/(2.*SIGMA*SIGMA));
	}
	//ic1dim.print("ic1DIM:");

	//Fourier Transform.
	cx_vec ic1d_k =  fft(ic1dim);

	//Grab the normalization constant from the zeroeth mode.
	complex<double> nrml = ic1d_k(0);
	//cout << "Normsl, Squared: " << nrml << "\t" << endl;
		//<< real(nrm) << "\t" << imag(nrm) << endl;
	//cout << "Just Sum: " << sum(ic1dim) << endl;

	return real(nrml);
}

void flucts_to_file(string dirNme, vec arr, uword ij)
{

        const char *fileNameR;
        ofstream fileStreamR;

	string sT = to_string(ij);
	std::string s0r = dirNme + "/" + "arrFluct_" + sT + ".dat";
	fileNameR = s0r.c_str();
	

        fileStreamR.open(fileNameR);
        if ( !fileStreamR.is_open() )
        {
                cout << "fileStream for flucts is broken" << endl;
                cout << endl;
                exit(1);
        }
        else
        {

		//print to file.
		for(uword j=0; j<arr.n_elem; ++j)
		{
                fileStreamR << j << "\t" << arr(j) << endl;
		}

                fileStreamR << endl;
                fileStreamR.close();
	}

	return;
}


void probs_to_file(string dirNme, int iT, cx_vec psi)
{

        const char *fileNameR;
        ofstream fileStreamR;

	double tiMe = initTime + iT*ds;
	string sT = to_string(tiMe);

	std::string s0r = dirNme + "/" + "prob_" + sT + ".dat";
	fileNameR = s0r.c_str();


        fileStreamR.open(fileNameR);
        if ( !fileStreamR.is_open() )
        {
                cout << "fileStream for probs is broken" << endl;
                cout << endl;
                exit(1);
        }
        else
        {

		cx_vec probs(NR, fill::zeros); 

		//assign a new psi_{r,s} matrix for ease, sum probs.
		for(uword rr=0; rr<psi.n_elem; ++rr)
                {
		probs(rr) = psi(rr)*conj(psi(rr));    
                }

		//print each contracted probability to file.
		for(uword j=0; j<probs.n_elem; ++j)
		{
                fileStreamR << j << "\t" << real(probs(j)) << endl;
		}

                fileStreamR << endl;
                fileStreamR.close();
	}

	return;
}

void wavefunc_to_file(string dirNme, int iT, cx_vec psiL)
{

        const char *fileNameX;
        ofstream fileStreamX;

	double tiMe = initTime + iT*ds;
	string sT = to_string(tiMe);

	std::string s0x = dirNme + "/" + "psi_" + sT + ".dat";
	fileNameX = s0x.c_str();

        //snprintf(fileNameR,s0r);
        //snprintf(fileNameS,s0s);


        fileStreamX.open(fileNameX);
        if ( !fileStreamX.is_open() )
        {
                cout << "fileStream for waveFunc is broken" << endl;
                cout << endl;
                exit(1);
        }
        else
        {

                for(uword k=0; k<psiL.n_elem; ++k)
                {
                uword ss = k % NR;
                uword rr = (k-ss)/NR;
                //cout << "mod " << ss << "\t" << rr << endl;
                //double x = -0.5*CELLSIZE+dx*k;
                fileStreamX << k << "\t" << ss << "\t" << rr << "\t"
                        << real(psiL(k)) << "\t" << imag(psiL(k)) << endl;
                }

                fileStreamX << endl;
                fileStreamX.close();
	}

	return;
}

std::string create_output_dir(double wN, double bigK, string typ)
{
	//This whole function is a bunch of terrible garbage. But it seems to
	//work. At least currently. The conversions are simply awful.
	//std::string str0 = std::to_string(w21v);
	//std::string str1 = std::to_string(wKin);


	std::stringstream stream0;
	stream0 << std::fixed << std::setprecision(2) << wN;
	std::string str0 = stream0.str();

	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << bigK;
	std::string str2 = stream.str();

	std::string str = "e" + str0 + typ + "_K" + str2;
	cout << str << endl;
	std::string strxy;

	//Create main outputs directory.
	struct stat st = {0};
	if (stat(str.c_str(), &st) == -1) 
	{
	strxy = str;
	cout << endl;
	cout << "Creating directory for output: " << strxy << endl;
	cout << endl;
    	mkdir(strxy.c_str(), 0700);
	//return str;
	}
	else if (stat(str.c_str(), &st) == 0) 
	{
	strxy = "e" + str0 + typ + "_K" + str2;
	cout << endl;
	cout << "Creating directory for output: " << str << endl;
	cout << endl;
    	mkdir(strxy.c_str(), 0700);
	//return str;
	}

	return strxy;

}

std::string create_subdir_gen(string str, string subname)
{
	std::string unc0 = subname;
	//std::string correl0 = std::to_string("correl");
	//std::string probs0 = std::to_string("probs");

	std::string subdir1 = str + "/" + unc0;

	//Create main outputs directory.
	struct stat st1 = {0};
	if (stat(subdir1.c_str(), &st1) == -1) 
	{
    	mkdir(subdir1.c_str(), 0700);
	}

	return subdir1;

}



std::string create_subdir(string str, string subname, uword ij)
{
	std::string unc0 = subname;
	//std::string correl0 = std::to_string("correl");
	//std::string probs0 = std::to_string("probs");

	std::string str0 = std::to_string(ij);
	std::string subdir1 = str + "/" + unc0 + str0;

	//Create main outputs directory.
	struct stat st1 = {0};
	if (stat(subdir1.c_str(), &st1) == -1) 
	{
    	mkdir(subdir1.c_str(), 0700);
	}

	return subdir1;

}






