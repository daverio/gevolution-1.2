#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP

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

class Global_defect
{
private:
  string runID_;
  string path_;
  
  Lattice lat_;
  Lattice klat_;
  Field<double> phi_defect_;
  Field<double> pi_defect_;

  Field<double> pi_defect_prev_;
  Field<double> rho_;
  Field<double> P_;
  
  Field<Imag> rho_k_;
  
  PlanFFT<Imag> planrho_;
  
 

public:
  Global_defect(){;}
  Global_defect(string settings_filename)
  {
    initialize(settings_filename);
  }
  ~Global_defect()
  {

  }
  unsigned long int random_seed();
  void initialize(string settings_filename);
  void loadSettings(string settings_filename);

  void generate_initCond();
  void create_directory();
  
  void evolve();
  double Friedmann_eq(double a);
  void next_cosmology(); //this fuction compute a and adot_over_a
  void field_leapfrog();
  double potential(Site & x);
  double potentialprime(Site & x, int comp);
  
  
  template<typename T>
  void averageField(Field<T> &f, string filename,int val);
  void averagephidefect();
  void averagerhodefect();
  
  double modsqphi(Site & x);
  double modsqrho(Site & x);
  
  void compute_rho_P_();
  void compute_pk_();
  
  bool output_now();
  void output();



};

////////////////////////////////////////////////////////////////
////    Initialize
//
////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//  This function initialises the value 
//
// Parameters:
//     dx_ = the distance between the lattice points
//     dt_ = the time step taken
//     lat_ = the real space lattice
//     klat_ = the fourier space lattice
//     phi_defect = the scalar field of the defect
//     pi_defect = the phi dot term used for the leapfrog algorithm
//     rho_ = the T00 term of the defect
//     rho_k_ = the T00 term of the defect in real space
//     P_ = the Tii term of the defect
//     rho_c0 = the current critical density
//     
///////////////////////////////////////////////////////////////

void Global_defect::initialize(string settings_filename)
{
  loadSettings(settings_filename);


  dx_ = physicalSize_/latticeSize_;
  dt_ = courantFactor_ * dx_;

  lat_.initialize(3,latticeSize_,1);
  klat_.initializeRealFFT(lat_,0);
  
  phi_defect_.initialize(lat_,nComponents_);
  phi_defect_.alloc();

  pi_defect_.initialize(lat_,nComponents_);
  pi_defect_.alloc();

  pi_defect_prev_.initialize(lat_,nComponents_);
  pi_defect_prev_.alloc();

  rho_.initialize(lat_);
  rho_.alloc();
  
  rho_k_.initialize(klat_);
  rho_k_.alloc();
  planrho_.initialize(&rho_,&rho_k_);
//  planrho_.alloc();  
  
  P_.initialize(lat_);
  P_.alloc();
  
  rho_c0 = 3*H0*H0/(8*M_PI*G);

  generate_initCond();
  
  COUT<<"Saving the files to:"<<path_<<endl<<endl;
  COUT<< "The initial conditions set are:"<<endl;
  COUT<< "start time is = "<<t_start_<<endl;
  COUT<< " End time is = "<<t_end_<<endl;
  COUT<< "Time interval is = "<<dt_<<endl;
  COUT<< "Lattice interval is = "<< dx_<<endl<<endl;
  COUT<< "The value set for univ is: "<<univ<<endl;
  COUT<< "Omega matter is = "<< omega_m<<endl;
  COUT<< "Omega radiation is = "<< omega_r<<endl;
  COUT<< "Omega lambda is = "<< omega_lambda<<endl;
  COUT<< "H0 is = "<< H0 <<endl;
  COUT<< "initial scale factor is = "<< a_i<<endl;
  COUT<< "The G is ="<<G<<endl;
  COUT<<"The current critical density is ="<<rho_c0<<endl;
  COUT<<"The val of pi is="<<M_PI<<endl<<endl;
}

/////////////////////////////////////////////////////////////
//  loadSettings
//
/////////////////////////////////////////////////////////////
// This function loads the settings and it also creates a 
// directory with the given runID to save the outputs in 
// and also saves the initial values to a separate 
// file: " settingsfile.txt"
//
////////////////////////////////////////////////////////////
// Parameters:
//    
//    settings_filename = the file in which the initial params are defined
/////////////////////////////////////////////////////////////

void Global_defect::loadSettings(string settings_filename)
{

  COUT<< "loadSettings: reading settings from : "<<settings_filename<<endl<<endl;


  SettingsFile setfile;
  setfile.open(settings_filename, SettingsFile::autoCreate);

  setfile.read("runID",runID_);
  setfile.read("path",path_);
  setfile.read("nComponents",nComponents_);
  setfile.read("latticeSize",latticeSize_);
  setfile.read("physicalSize",physicalSize_);
  setfile.read("courantFactor",courantFactor_);
  setfile.read("t_start",t_start_);
  setfile.read("t_end",t_end_);

  setfile.read("lambda",lambda_);
  setfile.read("eta2",eta2_);

  setfile.read("friction coefficient 1",Fric_term_1);
  setfile.read("friction coefficient 2",Fric_term_2);
  setfile.read("Dissipation time end",t_dis);
  setfile.read("Dissipation friction coefficient 1",dis_fric_term_1);
  setfile.read("Dissipation friction coefficient 2",dis_fric_term_2);
  
  setfile.read("univ",univ);
  setfile.read("omega_r",omega_r);
  setfile.read("omega_m",omega_m);
  setfile.read("omega_lambda",omega_lambda);
  setfile.read("H0",H0);
  setfile.read("a_i",a_i);
  setfile.read("G",G);
  setfile.read("bin_numbers",binnos);
  setfile.close();
  
  if(parallel.rank() == 0)
  { 
    if (mkdir(path_.c_str (), 0777) == -1) 
    cerr << "There is an error in creating directory! The directory already exists:  " << strerror(errno) << endl; 
    else
    cout << "Directory for" <<" "<<runID_<<" created"<<endl<<endl ;
  }
  
  ofstream settingsfile;
  settingsfile.open (path_  + runID_+"_settingsfile.txt",ios::trunc);
  settingsfile << "runID = " <<runID_ <<endl << "nComponents =" << nComponents_ <<endl << "Lattice Size ="<<latticeSize_ <<endl<<"Physical Size ="<<physicalSize_<<endl<<"Courant Factor ="<< courantFactor_ <<endl<<"t_start ="<<t_start_<<endl<<"t_end ="<<t_end_<<endl<<"lambda ="<<lambda_<<endl<<"eta2 ="<<eta2_<<endl<<"omega_r ="<<omega_r<<endl<<"omega_m ="<<omega_m<<endl<<"omega_r ="<<omega_r<<endl<<"omega_lambda ="<<omega_lambda<<endl<<"bin_numbers ="<<binnos<<endl;
  settingsfile.close();
}

/////////////////////////////////////////////////////////////
//  random_seed
//
/////////////////////////////////////////////////////////////
// This function is used to generate the initial seed to 
// generate the phi_defect and pi_defect 
//
////////////////////////////////////////////////////////////

unsigned long int Global_defect::random_seed()
  {
   struct timeval tv;
   gettimeofday(&tv,0);
   return(tv.tv_sec + tv.tv_usec);
  }

/////////////////////////////////////////////////////////////
//  generate_initCond
//
/////////////////////////////////////////////////////////////
// This function generates the initial conditions of phi_defect
// and pi_defect.In this the initial conditions are 
// taken to be a gaussian field
//
////////////////////////////////////////////////////////////

void Global_defect::generate_initCond()
{

  // TODO: understand how to change seed with gsl....
  Site x(lat_);

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  
  gsl_rng_default_seed = random_seed();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(x.first();x.test();x.next())
  {
    double phiNorm2 = 0;
    //theta = gsl_rng_uniform (r);
    //phi_defect_(x,0) = eta2*sin(theta);
    //phi_defect_(x,1) = eta2*cos(theta);
    for(int c=0;c<nComponents_;c++)
    {
      phi_defect_(x,c) = gsl_ran_gaussian (r,1);
      phiNorm2 += phi_defect_(x,c)*phi_defect_(x,c);
    }
    double ratio =  sqrt(eta2_/phiNorm2);
    for(int c=0;c<nComponents_;c++)
    {
      phi_defect_(x,c) *= ratio;
      pi_defect_(x,c) = 0;
    }
  }
  
  gsl_rng_free (r);

  phi_defect_.saveHDF5(path_ + runID_ + "_phi_defect_initCond.h5");
  pi_defect_.saveHDF5(path_  + runID_ + "_pi_defect_initCond.h5");

  COUT<< "Initial Condition generated"<<endl<<endl;

}

/////////////////////////////////////////////////////////////
//  evolve
//
/////////////////////////////////////////////////////////////
// This function is the main function which has steps to 
// evolve the fields using the leapfrog algorithm
//
////////////////////////////////////////////////////////////

void Global_defect::evolve()
{
    
    COUT<< "Starting main loop"<<endl;
    step = 0;
    t_ = t_start_;
    averagephidefect();
    averagerhodefect();
    
    ofstream phifile;
    phifile.open (path_  + runID_+"_average_phi.txt",ios::trunc);
    phifile.close();
    
    ofstream rhofile;
    rhofile.open (path_  + runID_+"_average_rho.txt",ios::trunc);
    rhofile.close();
    
//    ofstream omegafile;
//    omegafile.open (path_  + runID_+"_omega_.txt",ios::trunc);
//    omegafile.close();
  
  while(t_ <= t_end_)
  {
    t_ += dt_;
    next_cosmology();
    field_leapfrog();
    compute_rho_P_();
    averagephidefect();
    averagerhodefect();
    
    if(output_now())
    {
      output();
    }
    step++;

  }
}

/////////////////////////////////////////////////////////////
//  field_leapfrog
//
/////////////////////////////////////////////////////////////
// This function has the leapfrog algorithm implemented on 
// the equation of motion of the global defect.
//
////////////////////////////////////////////////////////////
// Parameters:
//    
//    fric_term_1 = the coefficient used to change the friction term set using the settings file
//    fric_term_2 = the coefficient used to change the friction term set using the settings file
/////////////////////////////////////////////////////////////

void Global_defect::field_leapfrog()
{
  Site x(lat_);

  for(x.first();x.test();x.next())
  {
    for(int c = 0;c<nComponents_;c++)
    {

      phi_defect_(x,c) += dt_ * pi_defect_(x,c);
    }
  }

  phi_defect_.updateHalo(); //update the value of phi in the halo
  
  double fric_term;
  
  if(dissipation)
  {
      if (t_<gdefects.diss_end)
      {
        fric_term = gdefects.friction_coeff;
      }
      else
      {
        fric_term = 1;
      }
  }
  else
  {
      fric_term = 1;
  }

  double c1 = (1.0 - dt_ * (fric_term*adot_overa_)) / (1.0 + dt_ * (fric_term*adot_overa_));
  double c2 = dt_ / (1.0 + dt_ * adot_overa_);
  double a2 = a_*a_;

  //cout<<dt_<<" "<<adot_overa_<<" "<<a_<<endl;
  //cout<<c1<<" "<<c2<<" "<<a2<<endl;

  // put what is in pi in pi_defect_prev
  //.... we then switch the data between pi and pi_defect_prev:
  double * temp = pi_defect_prev_.data_;
  pi_defect_prev_.data_ = pi_defect_.data_;
  pi_defect_.data_ = temp;

  for(x.first();x.test();x.next())
  {
    for(int c = 0;c<nComponents_;c++)
    {

      double lapPhi = -6.0 * phi_defect_(x,c) ;
      for(int i = 0 ; i<3 ; i++)lapPhi += phi_defect_(x+i,c) + phi_defect_(x-i,c);
      lapPhi /= dx_*dx_;
      pi_defect_(x,c) = c1 * pi_defect_prev_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c) );
    }
  }

}

/////////////////////////////////////////////////////////////
//  modsqphi
//
/////////////////////////////////////////////////////////////
// This function computes.the modulus of phi_defect square
//
////////////////////////////////////////////////////////////

double Global_defect::modsqphi(Site &x)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
  return pow(phiNorm2,0.5);
}

/////////////////////////////////////////////////////////////
//  averagephidefect
//
/////////////////////////////////////////////////////////////
// This function computes the <phi_defect> and writes it to a file
//
////////////////////////////////////////////////////////////

void Global_defect::averagephidefect()
{
  Site x(lat_);
  double phisum_ = 0;
  for(x.first();x.test();x.next())
  {
    phisum_ += modsqphi(x);
  }
  parallel.sum(phisum_);
  phiavg_ = phisum_/pow(latticeSize_,3);
  
  if(parallel.rank() == 0)
	{
    ofstream phifile;
    phifile.open (path_ + runID_+"_average_phi_defect.txt",std::ios_base::app);
    phifile << t_<<" "<<phiavg_<<endl;
    phifile.close();
  }
  //averageField(phi_defect_,"/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_phi.txt",0);
}

/////////////////////////////////////////////////////////////
//  averagerhodefect
//
/////////////////////////////////////////////////////////////
// This function computes the <rho_> and writes it to a file
//
////////////////////////////////////////////////////////////

void Global_defect::averagerhodefect()
{
  Site x(lat_);
  double rhosum_ = 0;
  for(x.first();x.test();x.next())
  {
    rhosum_ += rho_(x);
  }
  parallel.sum(rhosum_);
  rhoavg_ = rhosum_/pow(latticeSize_,3);
  if(parallel.rank() == 0)
  {
    ofstream rhofile;
    rhofile.open (path_ + runID_ + "_average_rho_defect.txt",std::ios_base::app);
    rhofile << t_<<" "<<rhoavg_<<endl;
    rhofile.close();
  }
}

////////////////////////////////////////////////////
//  potential
///////////////////////////////////////////////////
// This function computes the potential of the scalar field
//
//////////////////////////////////////////////////

double Global_defect::potential(Site & x)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
  return lambda_ * ( phiNorm2 - eta2_) * ( phiNorm2 - eta2_) / 2.0;
}

////////////////////////////////////////////////
//  potentialprime
//
////////////////////////////////////////////////
// This function computes the derivative of potential
//
///////////////////////////////////////////////

double Global_defect::potentialprime(Site & x, int comp)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
  return 2.0 * lambda_ * ( phiNorm2 - eta2_) *  phi_defect_(x,comp);
}

////////////////////////////////////////////////
//  Friedmann_eq
//
////////////////////////////////////////////////
// This function computes the friedmann equation 
//
///////////////////////////////////////////////
//  Parameters:
//
//    omega_lambda = density fraction of dark energy
//    omega_r = density fraction of radiation
//    omega_m = density fraction of matter
///////////////////////////////////////////////

double Global_defect::Friedmann_eq(double a)
{
  return H0*pow((omega_r + omega_m*a + omega_lambda*(a*a*a*a)),0.5);
}

void Global_defect::next_cosmology()
{
  if(univ== "matter")
  {
    a_ = t_ * t_ / (t_end_*t_end_);
    adot_overa_  = 2 / t_;
  }

  if(univ=="lcdm")
  {
    double k1;
    double k2;
    double k3;
    double k4;
    
    double rho_c;

    k1 = dt_*Friedmann_eq(a_i);
    k2 = dt_*Friedmann_eq(a_i + k1/2);
    k3 = dt_*Friedmann_eq(a_i + k2/2);
    k4 = dt_*Friedmann_eq(a_i + k3);
    a_i += k1/6 + k2/3 +k3/3 +k4/6;
    a_ = a_i;
    adot_overa_ = Friedmann_eq(a_)/a_;

    rho_c = 3*(Friedmann_eq(a_)*Friedmann_eq(a_)/(a_*a_*a_*a_))/(8*M_PI*G);
    
    double omegam=0;
    double omegar=0;
    double omegal=0;
    
    omegam = omega_m*rho_c0/ (a_*a_*a_*rho_c);
    omegar = omega_r*rho_c0/ (a_*a_*a_*a_*rho_c);
    omegal = 1 - omegam - omegar;
    
    parallel.sum(omegam);
    parallel.sum(omegar);
    parallel.sum(omegal);
    
    if(parallel.rank() == 0)
    {
      ofstream omegafile;
      omegafile.open (path_ + runID_+"_omega_.txt",std::ios_base::app);
      omegafile << a_<<" "<<omegam<<" "<<omegar<<" "<<omegal<<endl;
      omegafile.close();
    }
  }
}

bool Global_defect::output_now()
{
  return step%100==0?true:false;
}


////////////////////////////////////////////////
//  output
//
////////////////////////////////////////////////
// This function writws the output 
//
///////////////////////////////////////////////

void Global_defect::output()
{
  COUT<<"outputing field at t="<<t_<<endl;
  string filename_end= int2string(t_,99999)+".h5";
  

  phi_defect_.saveHDF5(path_ + runID_ + "_phi_defect_" + filename_end);
  pi_defect_.saveHDF5(path_ + runID_ + "_pi_defect_" + filename_end);

  //computeT();
  rho_.saveHDF5(path_ + runID_ + "_rho_" + filename_end);
  compute_pk_();
}

////////////////////////////////////////////////
//  compute_rho_P
//
////////////////////////////////////////////////
// This function computes the T00 and Tii 
//
///////////////////////////////////////////////

void Global_defect::compute_rho_P_()
{
  Site x(lat_);

  double a2 = a_*a_;
  rhoavg_ = 0;
  for(x.first();x.test();x.next())
  {
    double mpidot = 0;
    double temp;
    double gradPhi2 = 0;
    //double phinorm2 = 0;
    for(int c=0;c<nComponents_;c++)
    {
      temp = (pi_defect_prev_(x,c)+pi_defect_(x,c))/2.0;
      mpidot = temp*temp;
      //phinorm2 = phi_defect_(x,c)*phi_defect_(x,c);
      for(int i = 0;i<3;i++)
      {
        temp = ( phi_defect_(x+i,c) - phi_defect_(x-i,c) ) / 2.0 / dx_;
        gradPhi2 += temp*temp;
      }
    }  
    rho_(x) = mpidot / 2.0 / a2 + potential(x) + gradPhi2  / 2.0 / a2;
    P_(x) = mpidot / 2.0 /a2 - potential(x) - gradPhi2 / 6.0 / a2;
  }
}

////////////////////////////////////////////////
//  compute_pk
//
////////////////////////////////////////////////
// This function computes the powerspectrum 
//
///////////////////////////////////////////////
//  Parameters:
//
//    rho_k_ = the field used to compute the pk
//    binnos = contains the no of bins in k: taken from settings file
//    physicalSize_ = the box size: taken from settings file 
///////////////////////////////////////////////

void Global_defect::compute_pk_()
{
  planrho_.execute(FFT_FORWARD); 
  string filename_end= path_ + runID_ + "_powerspectrum" + int2string(t_,99999)+".txt";
  output_powerSpectrum(rho_k_,
                          filename_end,
                          binnos,
                          physicalSize_,
                          false,
                          false,
                          true,
                          false);
}

#endif



//template<typename T>

//void Global_defect::averageField(Field<T> &f, string filename,int val)
//{
//	Site x(f.lattice());
//	T sum = 0;
//	
//	if (val =0)
//	{
//	  for(x.first();x.test();x.next())
//    	{
//      	sum += modsqphi(x);
//    	}
//  }
//  else
//	{
//	  for(x.first();x.test();x.next())
//    	{
//      	sum += modsqrho(x);
//    	}
//  }
//  
//  parallel.sum(sum);
//  
//	T ave = sum/ (latticeSize_*latticeSize_*latticeSize_);
//	
//	if(parallel.rank() == 0)
//	{
//		ofstream file;
//		file.open(filename,std::ios_base::app);
//		file<< t_<<" "<<ave<<endl;
//		file.close();
//	}
//}
