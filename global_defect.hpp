#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/time.h>
#include<thread>
#include<chrono>
#include <fstream>
#include <bits/stdc++.h>

#include "LATfield2.hpp"
using namespace LATfield2;
#include "powerSpectra.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"

class GlobalDefect
{
private:

    double dx_;
    Lattice * lat_;
    Lattice * klat_;
    Field<double> phi_defect_;
    Field<double> pi_defect_;

    metadata * sim_;
    defects_metadata * defects_sim_;

public:

  void initialize(Lattice * lat, Lattice * klat_, metadata * sim, defects_metadata * defects_sim);
  void generate_init_cond();
  void update_phi()
  void update_pi(double dt, double a, double adot_overa);

};

void initialize(Lattice * lat, Lattice * klat_, double dx, metadata * sim, defects_metadata * defects_sim)
{
  dx_ = dx;
  lat_ = lat;
  /*
  
  */
}
void generate_init_cond()
{

}
void update_phi()
{
  Site x(lat_); // why it does not compile!
  for(x.first();x.test();x.next())
  {
    for(int c = 0; c < defects_sim_.nComponents;c++)
    {
      phi_defect_(x,c) += dt_ * pi_defect_(x,c);
    }
  }
  phi_defect_.updateHalo(); //update the value of phi in the halo
}
void update_pi(double dt, double a, double adot_overa)
{
  Site x(lat_); // why it does not compile!
  double c1 = (1.0 - dt_ * (adot_overa_)) / (1.0 + dt_ * (adot_overa_));
  double c2 = dt_ / (1.0 + dt_ * adot_overa_);
  double a2 = a_*a_;

  // put what is in pi in pi_defect_prev
  //.... we then switch the data between pi and pi_defect_prev:
  double * temp = pi_defect_prev_.data_;
  pi_defect_prev_.data_ = pi_defect_.data_;
  pi_defect_.data_ = temp;

  for(x.first();x.test();x.next())
  {
    for(int c = 0;c<gdefects.nComponents;c++)
    {

      double lapPhi = -6.0 * phi_defect_(x,c) ;
      for(int i = 0 ; i<3 ; i++)lapPhi += phi_defect_(x+i,c) + phi_defect_(x-i,c);
      lapPhi /= dx_*dx_;
      pi_defect_(x,c) = c1 * pi_defect_prev_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c,gdefects) );
    }
  }

}


#endif
