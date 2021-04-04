#ifndef LOCAL_DEFECT_HPP
#define LOCAL_DEFECT_HPP

#include "defect_base.hpp"



// #include <math.h>

class LocalDefect : public DefectBase
{
public:

//Generic lattice field evolution variables
  int  dim; //Number of spatial dimensions
  int  N;   //Number of lattice points in each direction

  Real dt;  //Time increment per timestep
  Real tStart;

  Real t;
  // int timeStep;    //Time step number
  // int timeStep2;   //Second time step count (used eg. for extra diffusive steps between timesteps)

  int phase;        //which phase of evolution is current?


//Fields
  Field<Imag> phi;
  Field<Imag> pi;
  Field<Real> theta;
  Field<Real> epsilon;

  //Simulation-specific parameters
    double lambda;       //quarter of dimensionless potential coupling in action (time dependent)
    Real lambdaEnd;    //quarter of dimensionless potential coupling in action (at tEnd)
    Real sigma;          //VEV of scalar product of phi with itself
    Real q;            //charge or rather gauge coupling (time dependent)
    Real qEnd;
    bool qCoreIndependent;      // allows for independent core growth

    Real lambdaCoreGrowthIndexB;  //Power of a that scalar core grows with after tCoreGrowth (natural value = -1)
    Real qCoreGrowthIndexB;  //Power of a that gauge core grows with after tCoreGrowth (natural value = -1)


    //read this in setting (TODO)

    Real tCoreGrowth;
    Real coreGrowthIndexB;


    //Simulation-specific entities
    Real a;           //Scalefactor at timestep
    Real aHalfPlus;   //Scalefactor at timestep plus half step
    Real aHalfMinus;  //Scalefactor at timestep minus half step
    Real qHalfPlus;  //Charge at timestep plus half step
    Real qHalfMinus;  //Charge at timestep minus half step

    Real lambdaCoreGrowth;  //Value of lambda at tCoreGrowth
    Real qCoreGrowth;       //Value of charge at tCoreGrowth
    Real aCoreGrowth;       //Scalefactor at tCoreGrowth



//CONSTRUCTOR
  // LocalDefect(double dtau);
  LocalDefect();

//INITIALIZATION
  void initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim);
  // void initialize();

//CORE FUNCTIONS
  void evolve(double dtau, double a);

  //MISCELANEOUS FUNCTIONS
    // void first(double dtau);           //Sets all parameters (eg. time) to start values
    void start();           //Initializes parameters (eg. tEnd=-1 => tEnd=N*dx/2)
    // void next();            //Updates parameters (eg. time) to current values
    void nextCosmology(double a, double aHalfPlus, double aHalfMinus);

};



//===============================
//CONSTRUCTOR====================
//===============================

// LocalDefect::LocalDefect(double dtau_)
LocalDefect::LocalDefect()
{
  dim = 3;
  // N = 32;

  phase = 2;

  // dt = dtau_;

  tStart = 0.0;

  lambda = 1.0;
  lambdaEnd = 1.0;
  sigma = 1.0;
  q = 1.0;
  qEnd = 1.0;
  qCoreIndependent = false;


  //(TODO: read the following in settings file instead)

  tCoreGrowth = 366;
  coreGrowthIndexB = -1.0;


}

// void LocalDefect::initialize()
void LocalDefect::initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim)
{

  dx_ = dx;
  lat_ = lat;
  klat_ = klat;
  defects_sim_ = defects_sim;
  sim_ = sim;

  phi.initialize(*lat_,1);
  pi.initialize(*lat_,1);
  theta.initialize(*lat_, dim);
  epsilon.initialize(*lat_, dim);
  phi.alloc();
  pi.alloc();
  theta.alloc();
  epsilon.alloc();

}


//==================================
//CORE FUNCTIONS====================
//==================================

void LocalDefect::evolve(double dt_, double a)
{

  double dt = dt_;
  double dx = *dx_;

  Site x(*lat_);
  int i;
  int j;

  //Setup required parameters
  Real  aadt_aadxdx = pow( a / aHalfPlus , Real(2) ) * dt / (dx*dx);
  Real  aadt2D_aadxdx = aadt_aadxdx * Real(2) * dim;
  Real  aaaaldt_2aa = 0.5 * pow( a, Real(4) ) * lambda * dt / pow(aHalfPlus, 2);
  Real  daa_aa = ( pow(aHalfPlus, 2) - pow(aHalfMinus, 2) ) / pow(aHalfPlus, 2);
  Real  dqq_qq = ( pow(qHalfPlus, 2) - pow(qHalfMinus, 2) ) / pow(qHalfMinus, 2);
  Real  ss = sigma*sigma;
  Real  dtqq_dxdxqq = dt * (qHalfPlus*qHalfPlus) / (dx*dx * q*q);
  Real  aaqq2dt = a*a * qHalfPlus*qHalfPlus * 2.0 * dt;
  Imag  deltaPi;  //Non-damping term change in pi(x)
  Real* deltaEpsilon = new Real[dim]; //Non-damping term change in epsilon(x,i)

  //Evolve primary fields
  for( x.first(); x.test(); x.next() )
  {
    phi(x) += dt * pi(x);
    for(i=0; i<dim; i++) { theta(x,i) += dt * epsilon(x,i); }
  }
  phi.updateHalo();
  theta.updateHalo();
  //Evolve time derivative fields
  for( x.first(); x.test(); x.next() )
  {
    deltaPi = - ( aadt2D_aadxdx + aaaaldt_2aa * ( phi(x).norm() - ss ) ) * phi(x);

    for(i=0; i<dim; i++)
    {
       deltaPi += aadt_aadxdx * ( phi(x+i) * expi(theta(x,i)) + phi(x-i) * expi(-theta(x-i,i)) );
      deltaEpsilon[i] = -aaqq2dt * ( phi(x+i) * expi(theta(x,i)) * phi(x).conj() ).imag();

      for(j=(i==0); j<dim; j++, j+=(j==i))
      {
	deltaEpsilon[i] -= dtqq_dxdxqq * ( sin( theta(x+i,j)   - theta(x,j)   - theta(x+j,i) + theta(x,i) )
					 - sin( theta(x+i-j,j) - theta(x-j,j) - theta(x,i)  + theta(x-j,i) ) );
      }

      if(phase!=1)
      {
	epsilon(x,i) *= (1.0 + dqq_qq);          //self-reliant - must be first epsilon update
	epsilon(x,i) += deltaEpsilon[i];
      }
      else
      {
    	// epsilon(x,i) = deltaEpsilon[i] / ( diffusiveFactor * dt ); //Cancel multiplication of dt and also reduce by a factor
      COUT << "difusive setup not implemented" <<endl;
      }
    }

    if(phase!=1)
    {
      pi(x) *= (1.0 - daa_aa);      //self-reliant - must be first pi update
      pi(x) += deltaPi;
    }
    else
    {
      // pi(x) = deltaPi / (diffusiveFactor * dt); //Cancel multiplication of dt and also reduce by a factor
      COUT << "difusive setup not implemented" <<endl;
    }
  }

  delete[] deltaEpsilon;


}


//===================================
//MISCELANEOUS FUNCTIONS=============
//===================================

// void LocalDefect::first()
// {
//   t = tStart-dt;
//   timeStep = 0;
//   // timeStep2 = 0; //extra diffusive steps between timesteps
//
//   this->start();
// }
//


void LocalDefect::start()
{

  // aCoreGrowth=bgc0.scaleFactor(tCoreGrowth);

  aCoreGrowth = 1. / (1. + defects_sim_->z_ic_defects);
  lambdaCoreGrowth = lambdaEnd * pow( aCoreGrowth, -2*(1+lambdaCoreGrowthIndexB) );


  if (!qCoreIndependent) {
      // lambdaCoreGrowthIndexA = coreGrowthIndexA;
      lambdaCoreGrowthIndexB = coreGrowthIndexB;
      // qCoreGrowthIndexA = coreGrowthIndexA;
      qCoreGrowthIndexB = coreGrowthIndexB;
  }

  qCoreGrowth = qEnd * pow( aCoreGrowth, -(1+qCoreGrowthIndexB) );

}

// void LocalDefect::next()
// {
//
//
//   if(phase>1)
//   {
//     timeStep++;
//     t += dt;
//     timeStep2=0;
//   }
//   else
//   {
//     // timeStep2++;
//     // if(timeStep2==diffusiveFactor) { timeStep++; timeStep2=0; }
//     // t += dt/diffusiveFactor;
//     COUT << "difusive setup not implemented" <<endl;
//
//   }
//
//   double a;
//   a = 1. / (1. + sim.z_in);
//
//   this->nextCosmology(a);
//
// }


void LocalDefect::nextCosmology(double a, double aHalfPlus, double aHalfMinus)
{

  // Real frac;
  // Real fracPrev;
  // if(t<tDiffusive)
  //   {
  //     if(phase<1) { phase=1; COUT<<"Starting diffusive phase at t="<<t<<endl; }
  //
  //     a=aDissipation;
  //     adot_a=Real(0);
  //     lambda=lambdaDissipation;
  //     q=qDissipation;
  //     aHalfPlus = a;
  //     aHalfMinus = a;
  //     qHalfPlus = q;
  //     qHalfMinus = q;
  //   }
  //   else if(t<tDissipation)
  //   {
  //     if(phase<2) { phase=2; COUT<<"Starting dissipative phase at t="<<t<<endl; }
  //
  //     //Set at-step values
  //     a=aDissipation;
  //     adot_a=Real(0);
  //     lambda=lambdaDissipation;
  //     q=qDissipation;
  //
  //     //Ramp charge cf. Vincent DPhil
  //     if(dissipationChargeRamp)
  //     {
  //       //Calculate frac of damping phase completed
  //       frac = ( (t-tStart) ) / (tDissipation-tStart);
  //       fracPrev = ( (t-dt-tStart) ) / (tDissipation-tStart);
  //
  //       //Ramp charge up during first half
  //       if(frac<0.5)
  //       {
  // 	//Ramp charge upto target value (gets half way 0.25 through damping phase)
  // 	q *= 0.5 * ( 1 + tanh( 8 * (frac-0.25) ) );
  //       }
  //
  //       //ReGauss at end of ramping
  //       if(frac>0.5 && fracPrev<0.5) { this->reGauss(); }
  //     }
  //
  //     //Fake off-step a values
  //     aHalfPlus = a * sqrt( 1 + 0.5*dissipation*dt );
  //     aHalfMinus = a * sqrt( 1 - 0.5*dissipation*dt );
  //
  //     //Fake off-step q values
  //     qHalfPlus = q * sqrt( 1 - 0.5*dissipation*dt );
  //     qHalfMinus = q * sqrt( 1 + 0.5*dissipation*dt );
  //   }
  //   else if(t<tCoreGrowth)
  //   {
    //   if(phase<3) { phase=3; COUT<<"Starting coreGrowth A phase at t="<<t<<endl; }
    //
    //   //Set at-step values
    //   a = aCoreGrowth * pow(t/tCoreGrowth, eraA);
    //   adot_a = eraA / t;
    //   lambda = lambdaCoreGrowth * pow( a/aCoreGrowth, -2*(1+lambdaCoreGrowthIndexA) );
    //   q = qCoreGrowth * pow( a/aCoreGrowth, -(1+qCoreGrowthIndexA) );
    //
    //   //Set off-step values
    //   aHalfPlus  = aCoreGrowth * pow( (t+dt/2)/tCoreGrowth, eraA);
    //   aHalfMinus = aCoreGrowth * pow( (t-dt/2)/tCoreGrowth, eraA);
    //   qHalfPlus  = qCoreGrowth * pow( aHalfPlus/aCoreGrowth, -(1+qCoreGrowthIndexA) );
    //   qHalfMinus = qCoreGrowth * pow( aHalfMinus/aCoreGrowth, -(1+qCoreGrowthIndexA) );
    // }
    // else
    // {
    //   if(phase<4){ phase=4; COUT<<"Starting coreGrowth B phase at t="<<t<<endl; }
    //   if (eraBcosmologyType == 1)
    //     {	    if(eraBparams[7]==0)
  	//     {
  	// 	  eraBparams[7]=1;
  	// 	  COUT<<"-------------------------------------"<<endl;
  	// 	  COUT<<"Starting Radiation-Matter transition."<<endl;
  	// 	  COUT<<"tSwitch   = "<< eraBparams[1] << endl;
  	// 	  COUT<<"tEquality = "<< eraBparams[0] << endl;
  	// 	  COUT<<"-------------------------------------"<<endl;
  	//     }
    //     }
    //   else if (eraBcosmologyType == 2)
    //     {	    if(eraBparams[7]==0)
  	//     {
  	// 	  eraBparams[7]=1;
  	// 	  COUT<<"-------------------------------------"<<endl;
  	// 	  COUT<<"LCDM cosmology."<<endl;
  	// 	  COUT<<"tSwitch   = "<< eraBparams[0] << endl;
  	// 	  COUT<<"Omega Radiation at tEnd = "<< eraBparams[1] << endl;
  	// 	  COUT<<"Omega Matter at tEnd = "<< eraBparams[2] << endl;
  	// 	  COUT<<"-------------------------------------"<<endl;
  	//     }
    //     }

      //Set at-step values
      // a = bgc0.scaleFactor(t);
      // adot_a = bgc0.hubble(t);
    // Real aHalfPlus;   //Scalefactor at timestep plus half step
    // Real aHalfMinus;  //Scalefactor at timestep minus half step

      lambda = lambdaCoreGrowth * pow( a/aCoreGrowth, -2*(1+lambdaCoreGrowthIndexB) );
      q=qCoreGrowth * pow( a/aCoreGrowth, -(1+qCoreGrowthIndexB) );

      //Set off-step values
      // aHalfPlus = bgc0.scaleFactor(t+dt/2.0);
      // aHalfMinus = bgc0.scaleFactor(t-dt/2.0);

      qHalfPlus  = qCoreGrowth * pow( aHalfPlus/aCoreGrowth, -(1+qCoreGrowthIndexB) );
      qHalfMinus = qCoreGrowth * pow( aHalfMinus/aCoreGrowth, -(1+qCoreGrowthIndexB) );

      // }

}





#endif
