#ifndef DEFECT_BASE_HPP
#define DEFECT_BASE_HPP


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include "metadata.hpp"
#include "parser.hpp"


class DefectBase
{
protected:

  Lattice * lat_;
  Lattice * klat_;
  double *dx_;
  metadata * sim_;
  defects_metadata * defects_sim_;
  
public:
  virtual void initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim) {};
  virtual void generate_init_cond() {}; 
  virtual void update_phi(double *dt) {};
  virtual void update_pi(double *dt, double *a, double *adot_overa) {};
  virtual void writedefectSnapshots(string h5filename,const int snapcount) {};
  
};

#endif
