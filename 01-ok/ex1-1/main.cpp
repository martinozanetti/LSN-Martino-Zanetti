/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex1-1

#include <iostream>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include "../../mylib/mylib.h"
#include <cmath>

using namespace std;

//int main (int argc, char *argv[]){

int main (){

   int M=100000;        // Total number of throws
   int N=100;           // Number of blocks = number of experiments
   int L=int(M/N);      // Number of throws in each block, use for M a multiple of N
   
   double sum_prog = 0; // variabili di appoggio per la media
   double var_prog = 0;
   double sig_prog = 0; // variabili di appoggio per la std (sigma)
   double varsig_prog = 0;

   // ========= Andamenti media e std ==============

   Random rnd;  
   rnd.SetSeed();

   ofstream WriteResult;
   WriteResult.open("results.out"); 

   // Ciclo dei blocchi
   for(int i=0; i<N; i++){
   double mean = 0; 
   double varsig = 0; 

   // Ciclo dei lanci in un blocco
      for(int j=0; j<L; j++){ 
         double r = rnd.Rannyu();
         mean += r;
         varsig += pow(r-0.5,2);
      }
      mean = mean/L;          // media
      sum_prog += mean;       // somma medie
      var_prog += mean*mean;  // somma medie quadrate

      varsig = varsig/L;      // media delle sigma
      sig_prog += varsig;     // somma media delle sigma
      varsig_prog += varsig*varsig; // somma quadratica media delle sigma

      // salvo su file: media e suo errore; sigma e suo errore
      if (WriteResult.is_open()){
         WriteResult << sum_prog/(i+1) << " " << error(sum_prog/(i+1),var_prog/(i+1),i) << " " //
                     << sig_prog/(i+1) << " " << error(sig_prog/(i+1),varsig_prog/(i+1),i) <<endl;
      } else {cerr << "PROBLEM: Unable to open results.out" << endl;}
   }
   WriteResult.close();
 

   // ========= CHI^2 test ==============

   ofstream WriteChi;
   WriteChi.open("chi.out"); 

   Random rnd2;
   rnd2.SetSeed(); 

   M = 100; // Numero di intervalli
   L = 1e4; // Numeri di lanci per intervallo

   // Ciclo dei blocchi/esperimenti
   for(int k = 0; k<M; k++){
      double chi = 0;
      double ni[M] = {};

   // Ciclo dei lanci in singolo esperimento
      for(int j=0; j<L; j++){
         double r = rnd2.Rannyu();

   // Riempio i bin degli intervalli 0 \ 0.01 \ 0.02 \ ... \ 0.99 \ 1.00  
   // >>>> l'algoritmo può essere migliorato: basta fare r*100 modulo 100 (o qualcosa di simile)    
         for(int i = 0; i<M; i++){
            if(r >= 1.0/M*i and r < 1.0/M*(i+1)){
               ni[i]++;
            }
         }
      }
   
   // Ora che ho riempito ni, posso calcolare il primo chi^2
      for(int i = 0; i<M; i++){
         chi += pow(ni[i]-L/M,2)/(L/M);
      }

      // Salvo il risultato 
      if (WriteChi.is_open()){
         WriteChi << chi << endl; //
      } else cerr << "PROBLEM: Unable to open chi.out" << endl;

   // Torno all'inizio del ciclo e ripeto per l'esperimento successivo
   }
   
   WriteChi.close();

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

