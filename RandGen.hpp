#ifndef RANDGEN_HPP
#define RANDGEN_HPP

#include <cstdlib>
#include <cmath>

//#include "parallel2d.hpp"

class RandGen
{
public:
  //CONSTRUCTOR
  RandGen();
  RandGen(long seed, int type = 0, int shuffleSize = 1);

  //CORE FUNCTIONS
  void   initialize(long seed, int type = 0, int shuffleSize = 1);
  double generate();
  double generateN();

  //SHUFFLE TABLE
  void shuffleInit(int shuffleSize);
  void shuffleSelect();
  void shuffleUpdate();

  //GENERATORS
  void   initType0(long seed);
  void   initType1(long seed);
  void   initType2(long seed);

  double genType0();
  double genType1();
  double genType2();

private:
  int type_;
  long state_;
  long state2_;

  long* shuffleTable_;
  int   shuffleSize_;
  int   shufflePosition_;

} randGen;


//CONSTRUCTOR===================
RandGen::RandGen()
{
  type_ = -1;
}

RandGen::RandGen(long seed, int type, int shuffleSize)
{
  this->initialize(seed, type, shuffleSize);
}

//CORE FUNCTIONS===============
void RandGen::initialize(long seed, int type, int shuffleSize)
{
  type_ = type;

  //If seed=0 then randomize based upon current time
  if(seed==0) { seed=time(NULL); }

  //Allow for different seeds on different processes
  seed = seed + 1234*parallel.rank();

  //Initialize generator
  if(type_==0)      { this->initType0(seed); }
  else if(type_==1) { this->initType1(seed); }
  else if(type_==2) { this->initType2(seed); }
  else
    {
      COUT<<"RandGen: Random generator type not recognized"<<endl;
      COUT<<"RandGen: Exiting..."<<endl;
      exit(EXIT_FAILURE);
    }

  //initialize shuffle table
  this->shuffleInit(shuffleSize);
}

double RandGen::generate()
{
  double random;

  //Select from shuffle table
  this->shuffleSelect();

  //Apply generator
  if(type_==0)      { random=this->genType0(); }
  else if(type_==1) { random=this->genType1(); }
  else if(type_==2) { random=this->genType2(); }

  //Update shuffle table
  this->shuffleUpdate();

  return random;
}

//==============================
//SHUFFLE TABLE=================
//==============================

//Bays-Durham suffle table
//
//Note does nothing with genType0

void RandGen::shuffleInit(int shuffleSize)
{
  if(type_==0 && shuffleSize>1)
    {
      COUT<<"RandGen: Random generator type 0 requested with shuffle table"<<endl;
      COUT<<"RandGen: Switching off shuffle table..."<<endl;
      shuffleSize=1;
    }

  shuffleSize_=shuffleSize;

  if(shuffleSize_>1)
    {
      shuffleSize_=shuffleSize;

      shuffleTable_=new long[shuffleSize_];

      //Fill shuffle table using genType1
      for(int i=0; i<shuffleSize; i++)
	{
	  this->genType1();
	  shuffleTable_[i]=state_;
	}
    }
}

void RandGen::shuffleSelect()
{
  if(shuffleSize_>1)
    {
      shufflePosition_ = state_ % shuffleSize_;
      state_ = shuffleTable_[shufflePosition_];
    }
}

void RandGen::shuffleUpdate()
{
  if(shuffleSize_>1)
    {
      shuffleTable_[shufflePosition_] = state_;
    }
}

//==============================
//TYPE 0 RANDOM NUMBER GENERATOR
//==============================

//C Standard Library generator

void RandGen::initType0(long seed)
{
  srand(seed);
}

double RandGen::genType0()
{
  return double(rand())/double(RAND_MAX);
}

//==============================
//TYPE 1 RANDOM NUMBER GENERATOR
//==============================

// Minimal random number generator of Park and Miller
//
// Algorithm:- integers generated as: I(j+1) = (a * I(j)) % m
//           - divided by m to give real
//           - mod performed using Schrage factorization
//             (a*z) % m = a*(z % q) - r*floor(z/q)    if >= 0
//                       =    ""     -     ""      + m otherwise
// Period: 2^31-2 = 2.1e9
//
// Known problems: very small numbers are followed by
// small numbers since they are only multiplied by 2e4.
// Other more subtle serial correlations

void RandGen::initType1(long seed)
{
  //Note: if I(j) = 0 gives I(j+1) = 0 so not allowed
  state_ = seed;
}

double RandGen::genType1()
{
  long m = 2147483647;
  long a = 16807;
  long q = 127773;
  long r = 2836;

  long flr = state_/q;

  state_ = a*( state_ - q * flr ) - r * flr;
  if(state_<0){ state_ += m; }

  return double(state_)/double(m);
}


//==============================
//TYPE 2 RANDOM NUMBER GENERATOR
//==============================

// L'Ecuyer long period generator using two Park and Miller
// generators and adding results together modulus one of
// their individual moduli.
//
// Period: 2.3e18
//
// Known problems: none other than its period.

void RandGen::initType2(long seed)
{
  state_=seed;
  state2_=long(genType1()*2147483398 + 1);   //Start second generator at a random point (but set by first)
}

double RandGen::genType2()
{
  long m1 = 2147483563;
  long m2 = 2147483399;
  long a1 = 40014;
  long a2 = 40692;
  long q1 = 53668;
  long q2 = 52774;
  long r1 = 12211;
  long r2 = 3791;

  long flr;
  long combined;

  flr = state_/q1;
  state_= a1*( state_ - q1 * flr ) - r1 * flr;
  if(state_<0){ state_ += m1; }

  flr = state2_/q2;
  state2_ = a2*( state2_ - q2 * flr ) - r2 * flr;
  if(state2_<0){ state2_ += m2; }

  combined = state_ - state2_;
  if(combined<0) { combined += m1-1; }

  return double(combined)/double(m1);
}

//===================================
//NORMALLY DISTRIBUTED RANDOM NUMBERS
//===================================

//Based on that from Numerical Recipes in C using the Box-Muller method
//
//Generates points in a "unit" square till it lies within a unit circle
//then applies a coordinate transformation to obtain a normally
//distributed random number with unit variance and zero mean

double RandGen::generateN()
{
  double x;
  double y;
  double rSqu;
  do
  {
    x = 2.0 * this->generate() - 1.0;
    y = 2.0 * this->generate() - 1.0;
    rSqu = x*x + y*y;
  }
  while(rSqu>=1.0 || rSqu==0.0);

  return sqrt(-2.0*log(rSqu)/rSqu) * x;
}

#endif
