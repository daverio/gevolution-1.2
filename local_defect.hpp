#ifndef LOCAL_DEFECT_HPP
#define LOCAL_DEFECT_HPP

#include "defect_base.hpp"

class LocalDefect:public DefectBase
{
public:

//Generic lattice field evolution variables
  int  dim; //Number of spatial dimensions
  int  N;   //Number of lattice points in each direction

//Fields
  Field<Imag> phi;
  Field<Imag> pi;
  Field<Real> theta;
  Field<Real> epsilon;

//CONSTRUCTOR
  LocalDefect();

//INITIALIZATION
  // void initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim));
  void initialize();

//CORE FUNCTIONS
  void evolve();

};



//===============================
//CONSTRUCTOR====================
//===============================

LocalDefect::LocalDefect()
{
  dim = 3;
  N = 32;
}

// void LocalDefect::initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim)
void LocalDefect::initialize()
{

  // dx_ = dx;
  // lat_ = lat;
  // klat_ = klat;
  // defects_sim_ = defects_sim;
  // sim_ = sim;

  // phi.initialize(lat_);
  // pi.initialize(lat_);
  // theta.initialize(lat_, defects_sim_->nComponents);
  // epsilon.initialize(lat_, defects_sim_->nComponents);
  // phi.alloc();
  // pi.alloc();
  // theta.alloc();
  // epsilon.alloc();

}


//==================================
//CORE FUNCTIONS====================
//==================================

void LocalDefect::evolve()
{
  // Site x(lattice);
  int i;
  int j;
}




#endif
