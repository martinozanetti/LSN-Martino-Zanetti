/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

void Random :: SetSeed(){
   int seed[4];   
   int p1, p2;
   // --- initialize random number generator
   ifstream Primes("../../random/Primes"); // put 2 primes in p1, p2
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../../random/seed.in"); // put 4 seeds in seed vector
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            this->SetRandom(seed,p1,p2); // initialization core
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   return;
}

void Random :: SetPrimesCouple(int n){

   int p1, p2;
   ifstream Primes("../../random/Primes");
   if (Primes.is_open()){
      for (int i=0; i<1+n;i++)      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   n3 = p1;
   n4 = p2;

   return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exponential(double lambda) {
   double y=Rannyu();
   double x = 0;
   if(lambda==0){
      cerr << "Exponential distribution needs non-zero decay rate parameter (lambda)" << endl;
   }else{
      x=-1/lambda * log(1-y);
   }
   return x;
}

double Random :: Cauchy(double mean, double gamma) {
   double y=Rannyu();
   double x = 0;
   if(gamma==0){
      cerr << "Cauchy distribution needs non-zero width parameter (lambda)" << endl;
   }else{
      x= gamma *tan(M_PI*(y-0.5)) + mean;
   }
   return x;
}

/// da correggere la dipendenza del return da m e q
double Random :: Linear(double m, double q, double min, double max) {
   double y = Rannyu();
   double x = 0;
   if((max-min)*(m*max*0.5+q)-1 > 1e-10 ){ // nota: nel confronto tra double bisogna confrontare con una tolleranza
      cerr << "Also linear distribution must be normalized to 1" << endl;
   }else{
      x = 1.-sqrt(1-y);
   }
   return x;
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
