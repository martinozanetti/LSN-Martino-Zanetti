/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex1-2

#include <iostream>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include <cmath>

using namespace std;

int main (){
   
   Random rnd;   
      rnd.SetSeed();
   
   int M=1e4;    // Numero totale di lanci
   int N[4] = {1,2,10,100}; // Numero di ripetizioni dell'esperimento
   

   ofstream stand("stand.out");
      if(!stand.is_open()) cerr << "PROBLEM: Unable to open stand.out" << endl;
   ofstream expon("expon.out");
      if(!expon.is_open()) cerr << "PROBLEM: Unable to open expon.out" << endl;
   ofstream loren("loren.out");
      if(!loren.is_open()) cerr << "PROBLEM: Unable to open loren.out" << endl;
   
   // ciclo sul numero di ripetizioni dell'esperimento
   for(int i=0; i<4; i++){

         // inizializzo le somme 1 volta sola
         double accu_std = 0;
         double accu_exp = 0;
         double accu_lor = 0;
         
   // ciclo dei lanci in ogni esperimento
      for(int j=0; j<M; j++){
         
         // azzero le somme 
         if(j>0){
            accu_std = 0;
            accu_exp = 0;
            accu_lor = 0;
         }

   // ciclo per il calcolo delle somme: 
   // salvo valori diretti se N=1, le medie N-esime per N>1
         for(int k=0; k<N[i]; k++){
            accu_std += rnd.Rannyu();
            accu_exp += rnd.Exponential(1.);
            accu_lor += rnd.Cauchy(0.,1.);
         }

      // divido le somme -> medie
         accu_std /= double(N[i]);
         accu_exp /= double(N[i]);
         accu_lor /= double(N[i]);
      
      // le registro nei file di output
         stand << accu_std << endl;
         expon << accu_exp << endl;
         loren << accu_lor << endl;
      /* 
         in questo modo in ciascun file ho 4 colonne
         incolonnate in serie, di 10.000 dati ciascuna,
         ciascun dato corrispondente a una S_N.
      */

      }
   }
   
   stand.close();
   expon.close();
   loren.close();

   return 0;
}


/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
