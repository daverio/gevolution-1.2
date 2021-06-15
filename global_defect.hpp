#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP

#include "defect_base.hpp"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include<thread>
#include<chrono>
#include <fstream>
#include <bits/stdc++.h>

#include "LATfield2.hpp"
using namespace LATfield2;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"
#include "powerSpectra.hpp"
//#include "read_text.hpp"

//#include <dir.h>
//#include<conio.h>


class GlobalDefect:public DefectBase
{
private:

/*
	double dx_;
	Lattice * lat_;
	Lattice * klat_;
	metadata * sim_;
	defects_metadata * defects_sim_;
*/
	Field<float> phi_defect_;
	Field<float> pi_defect_;
	Field<float> pi_defect_prev_;

	Field<double> T00_defect_;
	Field<double>Tii_defect_;

	Field<Imag> T00_defect_k_;
	PlanFFT<Imag> planT00_;

	double lambda;
	int firststep;

public:
	void initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim);
	void initialize(double *dx);
	void generate_init_cond(string h5filename, string filename_phi, string filename_pi);
	void update_phi(double *dt);
	void update_pi(double *dt, double *a, double *adot_overa);
	void writedefectSnapshots(string h5filename,const int snapcount);
	void defects_stat_output(); 
	void compute_Tuv_defect(double a);
	void write_Tuv_defect(string h5filename, const int snapcount);

	unsigned long int random_seed();
	double potential(Site & x);
	double potentialprime(Site & x, int comp);
	double modsqphi(Site & x);
	double averagephi();
	double averagerhodefect();
	void compute_defect_pk_(string h5filename, const int count);


};


void GlobalDefect::initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim)
{
	dx_ = dx;
	lat_ = lat;
	klat_ = klat;
	defects_sim_ = defects_sim;
	sim_ = sim;

	phi_defect_.initialize(*lat_, defects_sim_->nComponents);
	phi_defect_.alloc();

	pi_defect_.initialize(*lat_, defects_sim_->nComponents);
	pi_defect_.alloc();

	pi_defect_prev_.initialize(*lat_, defects_sim_->nComponents);
	pi_defect_prev_.alloc();

	T00_defect_.initialize(*lat_);
	T00_defect_.alloc();

	Tii_defect_.initialize(*lat_);
	Tii_defect_.alloc();

	Tuv_defect_.initialize(*lat_,4,4,symmetric);
	Tuv_defect_.alloc();
	
	T00_defect_k_.initialize(*klat_);
	T00_defect_k_.alloc();

	planT00_.initialize(&T00_defect_,&T00_defect_k_);

	firststep = 0;
}


unsigned long int GlobalDefect::random_seed()
{
	struct timeval tv;
	gettimeofday(&tv,0);
	return(tv.tv_sec + tv.tv_usec);
}


void GlobalDefect::generate_init_cond(string h5filename, string filename_phi = "", string filename_pi = "")
{
	if(filename_phi == "" && filename_pi == "")
	{
		Site x(phi_defect_.lattice());
		const gsl_rng_type * T;
		gsl_rng * r;

		gsl_rng_env_setup();

		gsl_rng_default_seed = random_seed();

		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		for(x.first();x.test();x.next())
		{
			double phiNorm2 = 0;
			for(int c = 0; c < defects_sim_->nComponents; c++)
			{
				phi_defect_(x,c) = gsl_ran_gaussian (r,1);
				phiNorm2 += phi_defect_(x,c)*phi_defect_(x,c);
			}
			double ratio =  sqrt(defects_sim_->eta2/phiNorm2);
			for(int c = 0; c < defects_sim_->nComponents; c++)
			{
				phi_defect_(x,c) *= ratio;
				pi_defect_(x,c) = 0;
			}
		}
		gsl_rng_free (r);
		COUT << "Initial values of phi and pi are loaded!" << endl;
	}
	else if(filename_phi != "" && filename_pi != "")
	{
		phi_defect_.loadHDF5(filename_phi);
		pi_defect_.loadHDF5(filename_pi);
		COUT << "loading data done!" << endl;
	}

#ifdef EXTERNAL_IO
	COUT << "Currently defect snapshot does not work with external IO" << endl;
#else
	phi_defect_.saveHDF5(h5filename + "phi_defect_init.h5");
	pi_defect_.saveHDF5(h5filename + "pi_defect_init.h5");
#endif

	COUT << COLORTEXT_MAGENTA << "Initial Condition for defect is generated" << COLORTEXT_RESET << endl << endl;
}


void GlobalDefect::update_phi(double *dt)
{ 
	double *dt_ = dt;
	Site x(phi_defect_.lattice());
	for(x.first();x.test();x.next())
	{
	for(int c = 0; c < defects_sim_->nComponents;c++)
	{
		phi_defect_(x,c) += *dt_ * pi_defect_(x,c);
	}
	}
	phi_defect_.updateHalo(); 
}


double GlobalDefect::potentialprime(Site & x, 
                                        int comp)
{
	double phiNorm2 = 0;
	for(int i =0; i < defects_sim_->nComponents; i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
	return 2.0 * lambda * ( phiNorm2 - defects_sim_->eta2) *  phi_defect_(x,comp);
}


void GlobalDefect::update_pi(double *dt,
                                double *a, 
                                double *adot_overa)
{
	double *dt_ = dt;
	double *a_ = a;
	double *adot_overa_ = adot_overa;
	double friction_coefficient = 1;
	lambda = defects_sim_->lambda0/ *a / *a;
	Site x(pi_defect_.lattice()); 

	if(defects_sim_->dissipation)
	{
	if(*a < 1 / (1 + defects_sim_->diss_end) )
	{
		friction_coefficient = defects_sim_->friction_coeff;
	}
	else
		friction_coefficient = 1;
	}

	double c1 = (1.0 - *dt_ * friction_coefficient * (*adot_overa_)) / (1.0 + *dt_ * friction_coefficient * (*adot_overa_));
	double c2 = *dt_ / (1.0 + *dt_ * (*adot_overa_));
	double a2 = *a_ * *a_;

	float * temp = pi_defect_prev_.data_;
	pi_defect_prev_.data_ = pi_defect_.data_;
	pi_defect_.data_ = temp;

	for(x.first();x.test();x.next())
	{
	for(int c = 0; c < defects_sim_->nComponents; c++)
	{
		double lapPhi = -6.0 * phi_defect_(x,c) ;
		for(int i = 0 ; i<3 ; i++)lapPhi += phi_defect_(x+i,c) + phi_defect_(x-i,c);
		lapPhi /=   *dx_ * *dx_ * sim_->boxsize * sim_->boxsize ;
		pi_defect_(x,c) = c1 * pi_defect_prev_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c) );
	}
	}
}


double GlobalDefect::potential(Site & x)
{
	double phiNorm2 = 0;
	for(int i =0; i < defects_sim_->nComponents; i++) phiNorm2 += phi_defect_(x,i) * phi_defect_(x,i);
	return lambda * ( phiNorm2 - defects_sim_->eta2) * ( phiNorm2 - defects_sim_->eta2) / 2.0;
}

void GlobalDefect::compute_Tuv_defect(double a)
{
	Site x(phi_defect_.lattice());

	double a2 = a * a;
	double temp;
	double gradPhi[3];
	for(x.first();x.test();x.next())
	{
		double mpidot = 0;
		double gradPhi2 = 0;
//		double gradPhi2_ = 0;

		for(int c=0; c < defects_sim_->nComponents;c++)
		{
		if(firststep == 0)
		{
			temp = pi_defect_(x,c);
			mpidot += temp*temp;
			
		}
		else
		{
			temp = (pi_defect_prev_(x,c)+pi_defect_(x,c))/2.0;
			mpidot += temp*temp;
		}
		
		for(int i = 0;i<3;i++)
		{
			gradPhi[i] = ( phi_defect_(x+i,c) - phi_defect_(x-i,c) ) / 2.0 / *dx_ / sim_->boxsize; 
			temp = ( phi_defect_(x+i,c) - phi_defect_(x-i,c) ) / 2.0 / *dx_ / sim_->boxsize;
//			gradPhi2_ += temp*temp;
		}
		}

		gradPhi2 = (gradPhi[0]*gradPhi[0])+(gradPhi[1]*gradPhi[1])+(gradPhi[2]*gradPhi[2]);

//		T00_defect_(x) = a2 * (mpidot / 2.0 / a2   + potential(x) + gradPhi2_ / 2.0 / a2) ;
		Tuv_defect_(x, 0, 0) = a2 * (mpidot / 2.0 / a2 + potential(x) + gradPhi2 / 2.0 / a2) ;

		Tuv_defect_(x, 1, 1) = a2 * (mpidot / 2.0 - potential(x)* a2 - gradPhi2 / 2.0 + gradPhi[0] * gradPhi[0]);
		Tuv_defect_(x, 2, 2) = a2 * (mpidot / 2.0 - potential(x)* a2 - gradPhi2 / 2.0 + gradPhi[1] * gradPhi[1]);
		Tuv_defect_(x, 3, 3) = a2 * (mpidot / 2.0 - potential(x)* a2 - gradPhi2 / 2.0 + gradPhi[2] * gradPhi[2]);

		Tuv_defect_(x, 1, 0) = a2 * gradPhi[0]*sqrt(mpidot);
		Tuv_defect_(x, 2, 0) = a2 * gradPhi[1]*sqrt(mpidot);
		Tuv_defect_(x, 3, 0) = a2 * gradPhi[2]*sqrt(mpidot);
		Tuv_defect_(x, 2, 1) = a2 * gradPhi[1]*gradPhi[0];
		Tuv_defect_(x, 3, 1) = a2 * gradPhi[2]*gradPhi[0];
		Tuv_defect_(x, 3, 2) = a2 * gradPhi[2]*gradPhi[1];

	}
	firststep = 1;
}

void GlobalDefect::write_Tuv_defect(string h5filename, const int snapcount)
{
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%05d", snapcount);

#ifdef EXTERNAL_IO
	COUT << "Currently defect snapshot does not work with external IO" << endl;
#else
	Tuv_defect_.saveHDF5(h5filename + filename_def + "_Tuv_global_defect_.h5" );
//	T00_defect_.saveHDF5(h5filename + filename_def + "_T00_global_defect_.h5" );
#endif
}


void GlobalDefect::writedefectSnapshots(string h5filename,
                                          const int snapcount)
{
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%03d", snapcount);

#ifdef EXTERNAL_IO
	COUT << "Currently defect snapshot does not work with external IO" << endl;
#else
	phi_defect_.saveHDF5(h5filename + filename_def + "_phi_defect_.h5");
	pi_defect_.saveHDF5(h5filename + filename_def + "_pi_defect_.h5");
#endif

}



double GlobalDefect::modsqphi(Site &x)
{
	double phiNorm2 = 0;
	for(int i = 0; i < defects_sim_->nComponents; i++ ) phiNorm2 += phi_defect_(x,i) * phi_defect_(x,i) ;
	return pow(phiNorm2,0.5);
}

double GlobalDefect::averagephi()
{
	Site x(phi_defect_.lattice());
	double phisum_ = 0;
	for(x.first();x.test();x.next())
	{
		phisum_ += modsqphi(x);
	}
	parallel.sum(phisum_);

	double size = sim_->numpts;
	double phiaverage_ =  phisum_/pow(size,3); 
	return phiaverage_; 
}


void GlobalDefect::defects_stat_output()
{
	double val = averagephi(); 
	COUT << " The average value of the field is = " << COLORTEXT_MAGENTA << val << COLORTEXT_RESET << endl; 
}

double GlobalDefect::averagerhodefect()
{
	Site x(phi_defect_.lattice());
	double rhoavg_;
	double rhosum_ = 0;
	double latsize = sim_->numpts;
	double lat3 = latsize*latsize*latsize;

	for(x.first();x.test();x.next())
	{
		rhosum_ += Tuv_defect_(x, 0, 0);
	}
	parallel.sum(rhosum_);
	rhoavg_ = rhosum_/lat3;
	return rhoavg_;
}

void GlobalDefect::compute_defect_pk_(string h5filename, const int count)
{
	planT00_.execute(FFT_FORWARD); 
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%03d", count);

	output_powerSpectrum(T00_defect_k_,
                          h5filename + filename_def + "defect_pk.txt",
                          sim_->numbins,
                          1.0,
                          false,
                          false,
                          true,
                          false);
}
#endif
