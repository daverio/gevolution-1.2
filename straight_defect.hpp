#ifndef STRAIGHT_DEFECT_HPP
#define STRAIGHT_DEFECT_HPP

#include "defect_base.hpp"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <thread>
#include <chrono>
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


class StraightDefect:public DefectBase
{
private:

	double v;
	double x_;
	double x0;
	double y0;
	double gamma;
	double mu;
	
	double Tuv_straight_[4][4];
	
public:
	void initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim);

	void update_phi(double *dt);

	void compute_Tuv_defect(double a);

	void write_Tuv_defect(string h5filename, const int snapcount);
	
	void projection_Tuu_straight_defect(double a);
	
	void projection_Ti0_straight_defect(double a);
	
	void symtensorProjection_comm(Field<Real> * Tij);
	
	void write_straight_defect_position(string h5filename, const int snapcount, double z, double dtau, double tau, double a, double Hconf);
};

void StraightDefect::initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim)
{	
	dx_ = dx;
	lat_ = lat;
	klat_ = klat;
	defects_sim_ = defects_sim;
	sim_ = sim;
	
	x0 = defects_sim_->string_x_pos;
	y0 = defects_sim_->string_y_pos;
	v = defects_sim_->string_velocity;

	mu = 75.3154 / sim_->boxsize / sim_->boxsize; 
	 
	gamma = 1/sqrt(1 - v*v);
	COUT << " mu = " << mu << endl; 
	COUT << " String velocity = " << v << endl;
	COUT << " gamma = " << gamma << endl; 
	
	x0 /= sim_->boxsize;
	y0 /= sim_->boxsize;
	x_ = x0;
	
	Tuv_defect_.initialize(*lat_,4,4,symmetric);
	Tuv_defect_.alloc();
	
	Tuv_straight_[0][0] = mu * gamma / *dx / *dx ;  //   
	
	COUT << " Tuv_straight_ [0][0] = " << Tuv_straight_[0][0] << endl;
	COUT << " DX = " << *dx << endl;
	Tuv_straight_[1][1] = mu * gamma * v * v / *dx / *dx;
	Tuv_straight_[2][2] = 0;
	Tuv_straight_[3][3] = - mu / *dx /gamma / *dx; 

	Tuv_straight_[1][0] = mu * gamma * v / *dx  / *dx;
	Tuv_straight_[2][0] = 0;
	Tuv_straight_[3][0] = 0;
	Tuv_straight_[2][1] = 0;
	Tuv_straight_[3][1] = 0;
	Tuv_straight_[3][2] = 0;
	
	COUT << "initialization of straight defect done!" << endl;
}

void StraightDefect::update_phi(double *dt)
{ 	
	x_ += v * *dt;
	
//	COUT << "dt = " << *dt << " and new position = " << x_ << endl;
	
	if(x_>1)
	{
#ifdef LATFIELD2_HPP
		COUT << COLORTEXT_RED << "The string has reached the edge of the box, so aborting further computation now!!" << COLORTEXT_RESET << endl;
		parallel.abortForce();
#endif
	}

	if (x_ > 1) x_ -= 1;

}

void StraightDefect::compute_Tuv_defect(double a)
{	
	Site xTuv(Tuv_defect_.lattice());
	
	/// Empting Halo
	
	long sizeLocalGross[3];
	long sizeLocal[3];
	int comp = 10;
	int halo = Tuv_defect_.lattice().halo();
	
	for(int i=0;i<3;i++)
	{
		sizeLocal[i]=Tuv_defect_.lattice().sizeLocal(i);
		sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
	}
    
	for(int k=0;k<sizeLocalGross[2];k++)
    {
        for(int j=0;j<sizeLocalGross[1];j++)
        {
        	for(int i=0;i<sizeLocalGross[0];i++)
        	{
            	for(int c=0;c<comp;c++)Tuv_defect_(setIndex(sizeLocalGross,i,j,k),c) = 0;
        	}
        }
    }

	projection_Tuu_straight_defect(a);
	
	projection_Ti0_straight_defect(a);
	
	symtensorProjection_comm(&Tuv_defect_);	

}

void StraightDefect::write_Tuv_defect(string h5filename, const int snapcount)
{
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%05d", snapcount);

#ifdef EXTERNAL_IO
	COUT << "Currently defect snapshot does not work with external IO" << endl;
#else
	COUT << " String x = "  << x_  << " and y position = " << y0 << endl;
	Tuv_defect_.saveHDF5(h5filename + filename_def + "_Tuv_straight_defect_.h5" );
#endif
}

void StraightDefect::write_straight_defect_position(string h5filename, const int snapcount, double z, double dtau, double tau, double a, double Hconf)
{
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%05d", snapcount);
	if(snapcount ==0)
	{
		ofstream straightdefectfile;
		if(parallel.rank() == 0)
		{
			straightdefectfile.open (h5filename + "straight_defect.txt",ios::trunc);
			straightdefectfile << "#z" << " " << "#dt" << " " << "#tau" << " " << "#a" << " " << "#Hconf" << " " << "#x" << " " << "#snapcount"<< endl;
			straightdefectfile.close();
			
			straightdefectfile.open (h5filename + "straight_defect.txt",std::ios_base::app);
			straightdefectfile << z << " " << dtau << " " << tau << " " <<  a << " " <<  Hconf << " " << x_ * sim_->boxsize << " " << snapcount << endl;
			straightdefectfile.close();
		}
	}
	else
	{
		ofstream straightdefectfile;
		if(parallel.rank() == 0)
		{
			straightdefectfile.open (h5filename + "straight_defect.txt",std::ios_base::app);
			straightdefectfile << z << " " << dtau << " " << tau << " " <<  a << " " <<  Hconf << " " << x_ * sim_->boxsize << " " << snapcount << endl;
			straightdefectfile.close();
		}
	}
}


void StraightDefect::projection_Tuu_straight_defect(double a)
{	
	Site xTuv(Tuv_defect_.lattice());
	
    double rescalPos[3];
    double rescalPosDown[3];
    
    double string_pos[2];
    double string_ref_coord[2];
    
    string_pos[0] = x_;
    string_pos[1] = y0; 
    
    for (int i=0; i<2; i++) string_ref_coord[i] = (string_pos[i]-fmod(string_pos[i],*dx_))/ *dx_;
    
	COUT << " grid coord = " << string_ref_coord[0] << " " << string_ref_coord[1] << endl;	
	
	for(int i=0;i<sim_->numpts;i++)
	{
 		
		if(xTuv.setCoord(string_ref_coord[0], string_ref_coord[1], i))
		{	
			 
 			for(int k =0;k<2;k++)
		    {
		        rescalPos[k] = (string_pos[k] - string_ref_coord[k] * *dx_) / *dx_;
		        rescalPosDown[k] = 1.0l - rescalPos[k];
		    }
				
			for (int m=0; m<4; m++)
			{	

				Tuv_defect_(xTuv,m,m)        = rescalPosDown[0] * rescalPosDown[1] * Tuv_straight_[m][m] / a / a ;
				Tuv_defect_(xTuv+1,m,m)      = rescalPosDown[0] * rescalPos[1] * Tuv_straight_[m][m]  / a / a ;
				Tuv_defect_(xTuv+0,m,m)      = rescalPos[0] * rescalPosDown[1] * Tuv_straight_[m][m]  / a / a ;
				Tuv_defect_(xTuv+0+1,m,m)    = rescalPos[0] * rescalPos[1] * Tuv_straight_[m][m]  / a / a ;
		
			}
		}
	}
}

void StraightDefect::projection_Ti0_straight_defect(double a)
{	
	Site xTuv(Tuv_defect_.lattice());
	
    double rescalPos[3];
    double rescalPosDown[3];
    
    double string_pos[2];
    double string_ref_coord[2];
    
    string_pos[0] = x_;
    string_pos[1] = y0; 
    
    // CIC + NGP hybrid projection (same as T0i projection in gevolution)
    
    for (int i=0; i<2; i++) string_ref_coord[i] = (string_pos[i]-fmod(string_pos[i],*dx_))/ *dx_;

	
	for(int i=0;i<sim_->numpts;i++)
	{
 		
		if(xTuv.setCoord(string_ref_coord[0], string_ref_coord[1], i))
		{	
			 
 			for(int k =0;k<2;k++)
		    {
		        rescalPos[k]= (string_pos[k] - string_ref_coord[k] * *dx_) / *dx_;
		        rescalPosDown[k]= 1.0l - rescalPos[k];
		    }
			
			for (int m=1; m<4; m++)
			{	

				Tuv_defect_(xTuv,m,0)        = rescalPos[1] * Tuv_straight_[m][0]  / a / a;
				Tuv_defect_(xTuv+1,m,0)      = rescalPosDown[1] * Tuv_straight_[m][0]  / a / a ;

			}
		}
	}
}

void StraightDefect::symtensorProjection_comm(Field<Real> * Tij)
{
	
    if(Tij->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = Tij->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=Tij->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;

    int comp=10;
    iref = sizeLocalGross[0]-halo;
    
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,halo,j,k),c) += (*Tij)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }

    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++) bufferSend[c+comp*(i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }

    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);
	
    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

    //send halo in direction Z

//    //pack data
//    iref=sizeLocalGross[2]-halo;
//    for(int j=0;j<(sizeLocalOne[1]-2);j++)
//    {
//        for(int i=0;i<imax;i++)
//        {
//            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
//        }
//    }

//    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);
//	
//	COUT << " iref = " << iref << endl;

//    //unpack data

//    for(int j=0;j<(sizeLocalOne[1]-2);j++)
//    {
//        for(int i=0;i<imax;i++)
//        {
//            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
//        }
//    }

    free(bufferRec);
    free(bufferSend);
}


#endif
