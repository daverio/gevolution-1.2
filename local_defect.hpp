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

  Real dx;  //Lattice spacing
  Real dt;  //Time increment per timestep
  Real tStart;

  Real tau;
  int timeStep;    //Time step number
  int timeStep2;   //Second time step count (used eg. for extra diffusive steps between timesteps)

  Real timeStepEnd;


  int phase;        //which phase of evolution is current?


  //Lattice
  // Lattice lattice; (changed to defect_base *lat_)

  //Energy-momentum tensor
  Lattice     emLattice;        //Zero halo!
  Field<Real> emTensor;

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

    Real lambdaCoreGrowthIndexA;  //Power of a that scalar core grows with after tCoreGrowth (natural value = -1)
    Real lambdaCoreGrowthIndexB;  //Power of a that scalar core grows with after tCoreGrowth (natural value = -1)
    Real qCoreGrowthIndexA;  //Power of a that gauge core grows with after tCoreGrowth (natural value = -1)
    Real qCoreGrowthIndexB;  //Power of a that gauge core grows with after tCoreGrowth (natural value = -1)


    //Simulation-specific parameters (read from file)

    double eraA;
    // double eraB;
    double tCoreGrowth;
    double coreGrowthIndexA;
    double coreGrowthIndexB;

    int outputEMconsEvery=0;


    //Simulation-specific entities
    Real a;           //Scalefactor at timestep
    Real aHalfPlus;   //Scalefactor at timestep plus half step
    Real aHalfMinus;  //Scalefactor at timestep minus half step
    Real qHalfPlus;  //Charge at timestep plus half step
    Real qHalfMinus;  //Charge at timestep minus half step

    Real lambdaCoreGrowth;  //Value of lambda at tCoreGrowth
    Real qCoreGrowth;       //Value of charge at tCoreGrowth
    Real aCoreGrowth;       //Scalefactor at tCoreGrowth

    cosmology *cosmo;
    double fourpiG;



//CONSTRUCTOR
  // LocalDefect(double dtau);
  LocalDefect();

//INITIALIZATION
  void initialize(Lattice * lat, Lattice * klat, double dt_, double *dx, metadata * sim, defects_metadata * defects_sim, cosmology *cosmo_, double fourpiG_);

  void generateIC_defects_test();

  // void initialize();

//CORE FUNCTIONS
  void evolve(double tau, double dtau, double a, int timeStep_);



  //MISCELANEOUS FUNCTIONS
    // void first(double dtau);           //Sets all parameters (eg. time) to start values
    void start();           //Initializes parameters (eg. tEnd=-1 => tEnd=N*dx/2)
    // void next();            //Updates parameters (eg. time) to current values

    void defects_stat_output();

    void defects_EmConservation_output();

    void emConservationInit(string preString, string postString);




  private:

    void nextCosmology(double a, double tau, double dt_);

    void emCalc();
    void emtCalcSite(Site x,Real * EMT);

    void emConservationCommon();
    void emConservationUpdatePre();
    void emConservationUpdate();
    void emConservationUpdatePost();


    //Energy-momentum tensor conservation

    int       emConservationStatus_;
    string    emConservationFileName_;
    ofstream  emConservationFile_;

    //Output fields============================

    bool timeToDoIt(int step, int number, int startStep, int endStep, string spacing)
    int roundNearest(double input)


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

  // tStart = 0.0;
  // lambda = 1.0;
  lambdaEnd = 1.0;
  // sigma = 1.0;
  q = 1.0;
  qEnd = 1.0;
  qCoreIndependent = false;

  timeStep = 0;
  timeStep2 = 0; //extra diffusive steps between timesteps


  emConservationStatus_=0;


}

// void LocalDefect::initialize()
void LocalDefect::initialize(Lattice * lat, Lattice * klat, double dt_, double *dx_, metadata * sim, defects_metadata * defects_sim, cosmology *cosmo_, double fourpiG_)
{



  dx = *dx_;
  lat_ = lat;
  klat_ = klat;
  defects_sim_ = defects_sim;
  sim_ = sim;

  cosmo = cosmo_;
  fourpiG = fourpiG_;

  phi.initialize(*lat_,1);
  pi.initialize(*lat_,1);
  theta.initialize(*lat_, dim);
  epsilon.initialize(*lat_, dim);
  phi.alloc();
  pi.alloc();
  theta.alloc();
  epsilon.alloc();


  tStart = defects_sim_->tStart;
  // tEnd = defects_sim_->tEnd;
  eraA = defects_sim_->eraA;
  // eraB = defects_sim_->eraB;
  tCoreGrowth = defects_sim_->tCoreGrowth;
  coreGrowthIndexA = defects_sim_->coreGrowthIndexA;
  coreGrowthIndexB = defects_sim_->coreGrowthIndexB;
  sigma = defects_sim_->sigma_loc;
  lambda = defects_sim_->lambda_loc;

  // timeStepEnd = roundDown( (defects_sim_->tEnd - defects_sim_->tStart) / dt_ );
  timeStepEnd = floor( (defects_sim_->tEnd - defects_sim_->tStart) / dt_ );

  if( defects_sim_->numberEMconsOutputs > 0 )
    {
      outputEMconsEvery = timeStepEnd / defects_sim_->numberEMconsOutputs;
      if( outputEMconsEvery < 3 ) { outputEMconsEvery = 3; }
    }



  // if( defects_sim_->num_local_defect_output>0 )
  //   {

      int box[3];

      #ifdef EMT_CG

          box[0] = sim->numpts / 2;
          box[1] = sim->numpts / 2;
          box[2] = sim->numpts / 2;

          emLattice.initialize(dim, box, 0);
          emTensor.initialize(emLattice, 4, 4, LATfield2::symmetric);
          emTensor.alloc();
      #else



      box[0] = sim->numpts;
      box[1] = sim->numpts;
      box[2] = sim->numpts;

          // emLattice.initialize(dim, N, 0);
          emLattice.initialize(dim, box, 0);
          emTensor.initialize(emLattice, 4, 4, LATfield2::symmetric);
          emTensor.alloc();
      #endif



      this->emConservationInit(sim_->output_path, ".dat");


    // }

}

void LocalDefect::generateIC_defects_test()
{
  Site x(*lat_);
  //Generate epsilon site-by-site
    for( x.first(); x.test(); x.next() )
    {
        for(int i=0; i<dim; i++)
        {
            epsilon(x,i) = Real(0);
            theta(x,i) = Real(0);
        }
        pi(x) = Imag(0,0);

        phi(x).real() =  1 * cos(x.coord(0));
        phi(x).imag() =  1 * sin(x.coord(0));
    }

    //Update halos
phi.updateHalo();
// pi.updateHalo(); //No halo update required as never referenced off-site
theta.updateHalo();
// epsilon.updateHalo(); //No halo update required as never referenced off-site


    // Site X(emLattice);
    //
    // for( X.first(); X.test(); X.next() )
    // {
    //
    //   for (int i=0;i<4;i++)
    //   {
    //     for (int j=0;j<=i;j++)
    //     {
    //       emTensor(X,i,j) = 1;
    //     }
    //
    //   }
    //
    // }



}


//==================================
//CORE FUNCTIONS====================
//==================================

void LocalDefect::evolve(double tau, double dt_, double a, int timeStep_)
{

  timeStep = timeStep_;

  this->nextCosmology(a,tau, dt_);

  double dt = dt_;
  // double dx = *dx_;

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


// COUT << defects_sim_->eraA << "is the era" << endl;

  // aCoreGrowth=bgc0.scaleFactor(tCoreGrowth);

  aCoreGrowth = 1. / (1. + defects_sim_->z_ic_defects);
  lambdaCoreGrowth = lambdaEnd * pow( aCoreGrowth, -2*(1+lambdaCoreGrowthIndexB) );


  if (!qCoreIndependent) {
      lambdaCoreGrowthIndexA = coreGrowthIndexA;
      lambdaCoreGrowthIndexB = coreGrowthIndexB;
      qCoreGrowthIndexA = coreGrowthIndexA;
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


void LocalDefect::nextCosmology(double a,double tau_, double dtau)
{

  tau = tau_;

  Real frac;
  Real fracPrev;
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
    // else if(t<tCoreGrowth)
    if(tau<tCoreGrowth)
    {


      if(phase<3) { phase=3; COUT<<"Starting coreGrowth A phase at t="<<tau<<endl; }

      //Set at-step values
      a = aCoreGrowth * pow(tau/tCoreGrowth, eraA);
      // adot_a = eraA / t;  //NOT NEEDED
      lambda = lambdaCoreGrowth * pow( a/aCoreGrowth, -2*(1+lambdaCoreGrowthIndexA) );
      q = qCoreGrowth * pow( a/aCoreGrowth, -(1+qCoreGrowthIndexA) );

      //Set off-step values
      aHalfPlus  = aCoreGrowth * pow( (tau+dtau/2)/tCoreGrowth, eraA);
      aHalfMinus = aCoreGrowth * pow( (tau-dtau/2)/tCoreGrowth, eraA);
      qHalfPlus  = qCoreGrowth * pow( aHalfPlus/aCoreGrowth, -(1+qCoreGrowthIndexA) );
      qHalfMinus = qCoreGrowth * pow( aHalfMinus/aCoreGrowth, -(1+qCoreGrowthIndexA) );
    }
    else
    {
      if(phase<4){ phase=4; COUT<<"Starting coreGrowth B phase at t="<<tau<<endl; }
      // if (eraBcosmologyType == 1)
      //   {	    if(eraBparams[7]==0)
  	  //   {
  		//   eraBparams[7]=1;
  		//   COUT<<"-------------------------------------"<<endl;
  		//   COUT<<"Starting Radiation-Matter transition."<<endl;
  		//   COUT<<"tSwitch   = "<< eraBparams[1] << endl;
  		//   COUT<<"tEquality = "<< eraBparams[0] << endl;
  		//   COUT<<"-------------------------------------"<<endl;
  	  //   }
      //   }
      // else if (eraBcosmologyType == 2)
      //   {	    if(eraBparams[7]==0)
  	  //   {
  		//   eraBparams[7]=1;
  		//   COUT<<"-------------------------------------"<<endl;
  		//   COUT<<"LCDM cosmology."<<endl;
  		//   COUT<<"tSwitch   = "<< eraBparams[0] << endl;
  		//   COUT<<"Omega Radiation at tEnd = "<< eraBparams[1] << endl;
  		//   COUT<<"Omega Matter at tEnd = "<< eraBparams[2] << endl;
  		//   COUT<<"-------------------------------------"<<endl;
  	  //   }
      //   }

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

      }

}



//===================================
//STAT FUNCTIONS=============
//===================================



void LocalDefect::defects_EmConservation_output()
{
  // double val = averagephi();
  // COUT << " The average value of the field is = " << COLORTEXT_MAGENTA << val << COLORTEXT_RESET << endl;


  // this->emCalc();

  //EM conservation testing==================

  // if(sim.timeStep+2<=timeStepEnd && (sim.timeStep+1)%outputEMconsEvery==0 && sim.timeStep2==0)
  if(timeStep+2<=timeStepEnd && (timeStep+1)%outputEMconsEvery==0 && timeStep2==0)
  {
    // if(timerEMcalc<Real(0)) { timerSplit=clock(); }

    // COUT<<"First EM conservation at               timestep: "<<sim.timeStep<<" time: "<<sim.t<<endl;
    COUT<<"First EM conservation at               timestep: "<<timeStep<<" time: "<<tau<<endl;
    this->emCalc();
    this->emConservationUpdatePre();

    // if(timerEMcalc<Real(0))
    // {
    //   timerEMcalc=seconds(clock(),timerSplit);
    //   runTimeEMcons = timerEMcalc * 3 * numberEMconsOutputs;
    //   COUT<<"Estimated EM calculation time: "<<hourMinSec(runTimeEMcons)<<endl;
    //   hibernateTime -= timerEMcalc/60; //Just in case it takes longer than grace period
    //   COUT<<"hibernateTime now " << hibernateTime << endl;
    // }
  // }
  // else if(sim.timeStep>0 && sim.timeStep+1<=timeStepEnd && sim.timeStep%outputEMconsEvery==0 && sim.timeStep2==0)
  }
  else if(timeStep>0 && timeStep+1<=timeStepEnd && (timeStep)%outputEMconsEvery==0 && timeStep2==0)
  {
    // COUT<<"Second EM conservation at              timestep: "<<sim.timeStep<<" time: "<<sim.t<<endl;
    COUT<<"Second EM conservation at              timestep: "<<timeStep<<" time: "<<tau<<endl;
    this->emCalc();
    this->emConservationUpdate();
  }
  // else if(sim.timeStep>1 && (sim.timeStep-1)%outputEMconsEvery==0 && sim.timeStep2==0)
  else if(timeStep>1 && (timeStep-1)%outputEMconsEvery==0 && timeStep2==0)
  {
    // COUT<<"Third EM conservation at               timestep: "<<sim.timeStep<<" time: "<<sim.t<<endl;
    COUT<<"Third EM conservation at               timestep: "<<timeStep<<" time: "<<tau<<endl;
    this->emCalc();
    this->emConservationUpdatePost();
  }






}



//ENERGY-MOMENTUM CONSERVATION OPERATIONS=====

void LocalDefect::emConservationInit(string preString, string postString)
{
  if( parallel.isRoot() )
    {
      emConservationFileName_=preString+"emConserve"+postString;

  //     if(hibernateWaking_==0)
	// {
	  emConservationFile_.open( emConservationFileName_.c_str(), fstream::out | fstream::trunc);
	// }
  //     else
	// {
	//   emConservationFile_.open( emConservationFileName_.c_str(), fstream::out | fstream::app);
	// }

      if(!emConservationFile_.good())
	{
	  cout<<"Simulation::emConvservationInit() - Could not open file: "<<emConservationFileName_<<endl;
	  cout<<"Simulation::emConvservationInit() - Exiting..."<<endl;
	  parallel.abortRequest();
	}

  //     if(hibernateWaking_==0)
	// {
	  emConservationFile_<<"#Energy-momentum conservation test"<<endl;
	  emConservationFile_<<"#Column 1: Time of Pre-step"<<endl;
	  emConservationFile_<<"#Column 2: Average EM-density at Pre-step"<<endl;
	  emConservationFile_<<"#Column 3: Time of Mid-step"<<endl;
	  emConservationFile_<<"#Column 4: adot-over-a at Mid-step"<<endl;
	  emConservationFile_<<"#Column 5: Average EM-density at mid-step"<<endl;
	  emConservationFile_<<"#Column 6: Average 3*pressure mid-step"<<endl;
	  emConservationFile_<<"#Column 7: Time of Post-step"<<endl;
	  emConservationFile_<<"#Column 8: Average EM-density at Post-step"<<endl;
	// }

      emConservationFile_.close();


      //Set status
      emConservationStatus_=1;
    }
  parallel.barrier();


}




#ifdef EMT_CG
void LocalDefect::emtCalcSite(Site x,Real * EMT)
{
    Real   T;
    Real   aaL;
    Imag*  D;
    Real** F;
    Real   EiEi;
    Real   FijFij;
    Real   DiDi;

    //Allocate memory
    //(for Fij, only j>i is allocated, so reference as F[i][j-i-1])
    D = new Imag[dim];
    F = new Real*[dim];
    for(int i=0; i<dim; i++) { F[i] = new Real[dim]; }

    //Calculate constants
    Real aaqdx = a*a * q * dx;
    Real aaqqdxdx = a*a * q*q * dx*dx;
    Real aa4 = a*a * Real(4);
    Real qdxdx = q * dx*dx;
    Real aal_4 = a*a * lambda / Real(4);
    Real ss = sigma*sigma;


    int diag[3];
    diag[0]=4;
    diag[1]=7;
    diag[2]=9;

    //cout << "Break point 1 passed by Process : "<< parallel.rank() <<endl;
    //parallel.barrier();


    //Calculate gauge derivative, field strength tensor, etc.
    EiEi = Real(0);
    FijFij = Real(0);
    DiDi = Real(0);
    for(int i=0; i<dim; i++)
    {
	EiEi += pow( epsilon(x,i), Real(2) );

	D[i] = ( phi(x+i)*expi(theta(x,i)) - phi(x) ) / dx;
	DiDi += D[i].norm();

	F[i][i] = Real(0);
	for(int j=i+1; j<dim; j++)
	{
	    F[i][j] = (Real(2)/qdxdx) * sin( (theta(x+i,j) - theta(x,j) - theta(x+j,i) + theta(x,i)) / Real(2) );
	    F[j][i] = -F[i][j];
	    FijFij += pow( F[i][j], Real(2) );
	}
    }
    FijFij *= Real(2); //Double because sum above is only over half of terms

    //cout << "FijFij passed by Process : "<< parallel.rank() <<endl;
    //parallel.barrier();


    //Calculate Legrangian density
    aaL = 0.5*EiEi/aaqqdxdx - FijFij/aa4 + pi(x).norm() - DiDi;
    aaL -= aal_4*pow( phi(x).norm() - ss, Real(2) );

    //CALCULATE EM TENSOR COMPONENTS (doubly lowered indices)

    //Time-Time component
    EMT[0] = EiEi/aaqqdxdx + Real(2) * pi(x).norm() - aaL;

//Time-Space components
    for(int i=0; i<dim; i++)
    {
	T = Real(0);
	for(int j=0; j<dim; j++)
	{
	    T += epsilon(x,j) * F[i][j];
	}
	EMT[1+i] = T / aaqdx + Real(2) * ( pi(x).conj()*D[i] ).real();
    }

    //Space-Space components (off-diagonal)
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<i; j++)
	{
            T = Real(0);
            for(int k=0; k<dim; k++) { T += F[i][k] * F[j][k]; }
	    int index;
	    index = int( abs(i-j) + (j+1)*(4+0.5-0.5*(j+1) ) );

	    EMT[index] = -epsilon(x,i)*epsilon(x,j)/aaqqdxdx + T/(a*a) + Real(2)*(D[i].conj()*D[j]).real();
	}
    }
    //Space-Space components (diagonal)
    for(int i=0; i<dim; i++)
    {
	T = Real(0);
	for(int k=0; k<dim; k++) { T += pow( F[i][k], Real(2) ); }

	EMT[diag[i]] = -pow(epsilon(x,i),Real(2))/aaqqdxdx + T/(a*a) + Real(2)*D[i].norm() + aaL;
    }


    delete D;
    delete F;


}
#endif

void LocalDefect::emCalc()
{

  #ifdef EMT_CG

      Real T[10];
      Real em[10];

      for(int i =0 ;i<10;i++)T[i]=0;

      Site x(*lat_);
      Site X(emLattice);


  //    x.first();
  //    this->emtCalcSite(x,T);

  //    for(int i =0 ;i<10;i++)COUT<<"T["<<i<<"] = "<<T[i]<<endl;

      for(X.first();X.test();X.next())
      {

  	for(int i =0 ;i<10;i++)T[i]=0;

  	x.setCoord(2*X.coord(0),2*X.coord(1),2*X.coord(2));


  	this->emtCalcSite(x,T);
  	for(int i =0 ;i<10;i++)em[i]=T[i];

  	this->emtCalcSite(x+0,T);
          for(int i =0 ;i<10;i++)em[i]+=T[i];

  	this->emtCalcSite(x+1,T);
  	for(int i =0 ;i<10;i++)em[i]+=T[i];

  	this->emtCalcSite(x+2,T);
  	for(int i =0 ;i<10;i++)em[i]+=T[i];

  	this->emtCalcSite(x+0+1,T);
  	for(int i =0 ;i<10;i++)em[i]+=T[i];

  	this->emtCalcSite(x+0+2,T);
  	for(int i =0 ;i<10;i++)em[i]+=T[i];

  	this->emtCalcSite(x+1+2,T);
  	for(int i =0 ;i<10;i++)em[i]+=T[i];

  	this->emtCalcSite(x+0+1+2,T);
  	for(int i =0 ;i<10;i++)em[i]+=T[i];


  	for(int i =0 ;i<10;i++)emTensor(X,i)=em[i]/=Real(8);
      }




  #else

  Site x(*lat_);
  Site X(emLattice);

  Real   T;
  Real   aaL;
  Imag*  D;
  Real** F;
  Real   EiEi;
  Real   FijFij;
  Real   DiDi;

  //Allocate memory
  //(for Fij, only j>i is allocated, so reference as F[i][j-i-1])
  D = new Imag[dim];
  F = new Real*[dim];
  for(int i=0; i<dim; i++) { F[i] = new Real[dim]; }

  //Calculate constants
  Real aaqdx = a*a * q * dx;
  Real aaqqdxdx = a*a * q*q * dx*dx;
  Real aa4 = a*a * Real(4);
  Real qdxdx = q * dx*dx;
  Real aal_4 = a*a * lambda / Real(4);
  Real ss = sigma*sigma;

  // cout << "Break point 1 passed by Process : "<< parallel.rank() <<endl;
  // parallel.barrier();



  for(x.first(), X.first(); x.test(); x.next(), X.next())
  {
    //Calculate gauge derivative, field strength tensor, etc.
    EiEi = Real(0);
    FijFij = Real(0);
    DiDi = Real(0);
    for(int i=0; i<dim; i++)
    {
      EiEi += pow( epsilon(x,i), Real(2) );

      D[i] = ( phi(x+i)*expi(theta(x,i)) - phi(x) ) / dx;
      DiDi += D[i].norm();

      F[i][i] = Real(0);
      for(int j=i+1; j<dim; j++)
      {
        F[i][j] = (Real(2)/qdxdx) * sin( (theta(x+i,j) - theta(x,j) - theta(x+j,i) + theta(x,i)) / Real(2) );
        F[j][i] = -F[i][j];
        FijFij += pow( F[i][j], Real(2) );
      }
    }
    FijFij *= Real(2); //Double because sum above is only over half of terms

    //cout << "FijFij passed by Process : "<< parallel.rank() <<endl;
    //parallel.barrier();


    //Calculate Legrangian density
    aaL = 0.5*EiEi/aaqqdxdx - FijFij/aa4 + pi(x).norm() - DiDi;
    aaL -= aal_4*pow( phi(x).norm() - ss, Real(2) );

    //CALCULATE EM TENSOR COMPONENTS (doubly lowered indices)

    //Time-Time component
    emTensor(X,0,0) = EiEi/aaqqdxdx + Real(2) * pi(x).norm() - aaL;

    //cout << "T-T component passed by Process : "<< parallel.rank() <<endl;
    //parallel.barrier();

    //Time-Space components
    for(int i=0; i<dim; i++)
    {
      T = Real(0);
      for(int j=0; j<dim; j++)
      {
        T += epsilon(x,j) * F[i][j];
      }
      emTensor(X,0,1+i) = T / aaqdx + Real(2) * ( pi(x).conj()*D[i] ).real();
    }

    //Space-Space components (off-diagonal)
    for(int i=0; i<dim; i++)
    for(int j=0; j<i; j++)
    {
      T = Real(0);
      for(int k=0; k<dim; k++) { T += F[i][k] * F[j][k]; }

      emTensor(X,1+i,1+j) = -epsilon(x,i)*epsilon(x,j)/aaqqdxdx + T/(a*a) + Real(2)*(D[i].conj()*D[j]).real();
    }

    //Space-Space components (diagonal)
    for(int i=0; i<dim; i++)
    {
      T = Real(0);
      for(int k=0; k<dim; k++) { T += pow( F[i][k], Real(2) ); }

      emTensor(X,1+i,1+i) = -pow(epsilon(x,i),Real(2))/aaqqdxdx + T/(a*a) + Real(2)*D[i].norm() + aaL;
    }
  }

  delete D;
  delete F;


  #endif


  // emOndate = true;
  // parallel.barrier();
}


void LocalDefect::emConservationUpdatePre()
{
  if(parallel.isRoot())
    {
      //Check status
      if(emConservationStatus_!=1)
	{
	  COUT<<"Simulation::emConvservationUpdatePre() - object not is ready state for this function"<<endl;
	  COUT<<"Simulation::emConvservationUpdatePre() - Exiting..."<<endl;
	  parallel.abortRequest();
	}

      this->emConservationCommon();

      //Write Pre time
      emConservationFile_<<tau<<" ";
    }
  parallel.barrier();

  //Sum current energy-density
  Site x(emLattice);
  Real T00Sum=0.0;

  for(x.first(); x.test(); x.next())
    {
      T00Sum += emTensor(x,0,0);
    }
  parallel.sum(T00Sum);

  if(parallel.isRoot())
    {
      //Write average energy density
      // emConservationFile_<<T00Sum/pow(N, Real(dim))<<" ";
      emConservationFile_<<T00Sum/pow(sim_->numpts, Real(dim))<<" ";

      emConservationFile_.close();
    }

  //Set status
  emConservationStatus_=2;
}

void LocalDefect::emConservationUpdate()
{
  //Check status
  if(emConservationStatus_!=2)
    {
      COUT<<"Simulation::emConvservationUpdate() - object not is ready state for this function"<<endl;
      COUT<<"Simulation::emConvservationUpdate() - Exiting..."<<endl;
      exit(EXIT_FAILURE);
    }

  if(parallel.isRoot())
    {
      this->emConservationCommon();

      //Write time and adot_a
      // emConservationFile_<<tau<<" "<<adot_a<<" ";
      emConservationFile_<<tau<<" "<<Hconf(a, fourpiG, *cosmo)<<" ";


    }

  //Sum current energy-density and pressure
  Site x(emLattice);
  Real T00Sum=0.0;
  Real TjjSum=0.0;

  for(x.first(); x.test(); x.next())
    {
      T00Sum += emTensor(x,0,0);
      for(int j=0; j<dim; j++) { TjjSum += emTensor(x,j+1,j+1); }
    }
  parallel.sum(T00Sum);
  parallel.sum(TjjSum);

  if(parallel.isRoot())
    {
      //Write averages
      // emConservationFile_<<T00Sum/pow(N, Real(dim))<<" "<<TjjSum/pow(N, Real(dim))<<" ";
      emConservationFile_<<T00Sum/pow(sim_->numpts, Real(dim))<<" "<<TjjSum/pow(sim_->numpts, Real(dim))<<" ";

      emConservationFile_.close();
    }

  //Set status
  emConservationStatus_=3;
}

void LocalDefect::emConservationUpdatePost()
{
  //Check status
  if(emConservationStatus_!=3)
    {
      COUT<<"Simulation::emConvservationUpdatePost() - object not is ready state for this function"<<endl;
      COUT<<"Simulation::emConvservationUpdatePost() - Exiting..."<<endl;
      exit(EXIT_FAILURE);
    }

   if(parallel.isRoot())
    {
      this->emConservationCommon();

      //Write Post time
      emConservationFile_<<tau<<" ";
    }

   //Sum current energy-density
  Site x(emLattice);
  Real T00Sum=0.0;

  for(x.first(); x.test(); x.next())
    {
      T00Sum += emTensor(x,0,0);
    }
  parallel.sum(T00Sum);

  if(parallel.isRoot())
    {
      //Write average energy density
      // emConservationFile_<<T00Sum/pow(N, Real(dim))<<endl;
      emConservationFile_<<T00Sum/pow(sim_->numpts, Real(dim))<<endl;

      emConservationFile_.close();
    }

  //Set status
  emConservationStatus_=1;
}


void LocalDefect::emConservationCommon()
{
  emConservationFile_.open( emConservationFileName_.c_str(), fstream::out | fstream::app);

  if(!emConservationFile_.good())
    {
      cout<<"Simulation::emConvservationCommon() - Could not open file: "<<emConservationFileName_<<endl;
      cout<<"Simulation::emConvservationCommon() - Exiting..."<<endl;
      exit(EXIT_FAILURE);
    }

  emConservationFile_.precision(7);
  emConservationFile_.width(emConservationFile_.precision()+7);
  emConservationFile_.setf(fstream::scientific | fstream::showpos);
}


//Output fields============================


void defects_stat_output()
{

  // if( timeToDoIt(sim.timeStep,numberFieldOutputs,timeStepFields,timeStepEnd,spacingFieldOutputs) && sim.timeStep2==0) // Old version saved fields 1 timestep different. MBH
  // 	{
  //
  //   }

}

bool LocalDefect::timeToDoIt(int step, int number, int startStep, int endStep, string spacing)
{
  double dlog;     //Increment in log(step/startStep) required
  int acceptStep;
  int lastAcceptStep = startStep;

  //If step is the startStep, then definitely accept
  if(step==startStep) { return 1; }

  if (spacing == "log"){

	if(startStep==0) { startStep = 1; } //Otherwise messes up following division

	//Calculate preliminary dlog value
	dlog = log10( double(endStep)/double(startStep) ) / double(number-1);

	for(int i=2; i<=number; i++)
	  {
		//Calculate nearest integer to log required
		acceptStep = roundNearest( lastAcceptStep * pow(double(10), dlog) );

		//If step unchanged, then increase by one
		if(acceptStep == lastAcceptStep) { acceptStep++; }

		//Store last step
		lastAcceptStep = acceptStep;

		//Redistribute logs over remaining steps
		dlog = log10( double(endStep)/double(acceptStep) ) / double(number-i);

		if(number > 0 && step > startStep && step == acceptStep)
		{
  //		COUT << " log";
		  return 1;
		}
	  }
  }
  else 	{
	  int outputEvery;

	  if (spacing != "linear"){
		COUT << "timeToDoIt: warning: spacing not log or linear ... linear assumed." << endl;
	  }

	  outputEvery = (endStep - startStep + 1) / number;

	  if (number > 0 && step > startStep && (step-startStep)%outputEvery==0 )
	  {
		//		COUT << " linear";
		return 1;
	  }
  }
  return 0;
}


int LocalDefect::roundNearest(double input)
{
  double f=floor(input);
  if(input-f<0.5) { return int(f); }
  else { return int(f)+1; }
}






#endif
