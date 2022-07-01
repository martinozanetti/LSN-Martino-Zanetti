/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex1-3

#include <iostream>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include "../../mylib/mylib.h"
#include <cmath>
#include <stdlib.h>
#include <iomanip>

using namespace std;

int main (){
   
   int M=1000;  // number of throws in each block
   int N=100;    // Number of blocks = number of experiments
   double L = 3; // needle Length
   double d = 5.; // Distance between lines
   double pi = 0.;

   double sum_prog = 0;
   double var_prog = 0;
   double N_hit = 0;
   double err = 0; 
 
   Random rnd;    
      rnd.SetSeed();
   
   ofstream Result("results.out");
      if(!Result.is_open()) cerr << "PROBLEM: Unable to open results.out" << endl;
   ofstream Result2("results2.out");
      if(!Result2.is_open()) cerr << "PROBLEM: Unable to open results2.out" << endl;
   

   // ciclo per l'aumento del numero di blocchi
   for(int i=0; i<N; i++){
      N_hit = 0;

   // ciclo degli esperimenti con M lanci ciascuno 
      for(int j=0; j<M; j++){

         // POSIZIONE DELLA CRUNA: y_pos casuale in [0,d) (uso p.b.c.!)
         double y_pos = d*rnd.Rannyu();

         // ANGOLO: x,y casuali in [0,+L)
         double x = L*rnd.Rannyu();
         double y = L*rnd.Rannyu();

         // normalizzazione: tengo solo i lanci nel cerchio, poi li porto sulla circonferenza
         if(y>sqrt(L*L-x*x)){ // lancio esterno invalido: torno allo step precedente
            j--;
         } 
         else if(y<=sqrt(L*L-x*x)){ // lancio interno valido: salvo i risultati
            // normalizzo ...  
            double a = x;
            x = x*L/sqrt(x*x+y*y);
            y = y*L/sqrt(a*a+y*y);

            // ... e salvo la posizione
            // (per vedere graficamente se è su una circonferenza)
            if (Result2.is_open()){
               Result2 << x << ' ' << y << endl;
            } else cerr << "PROBLEM: Unable to open results2.out" << endl;

            // se (y_pos+y >= L) o (<= 0) aggiungo un conteggio a N_hit
            if(y_pos+y >= d or y_pos+y <= 0) N_hit ++;
         }
      }

      // calcolo PI secondo Buffon
      pi = 2.0*L*double(M/N_hit)/d; 
      sum_prog += pi;     
      var_prog += pi*pi; 

      // salvo PI con il rispettivo errore
      // >>>> l'algoritmo può essere migliorato, basta prendere la parte intera di log10(err)
      if (Result.is_open()){
         err = 0;
         if (i>0) err = error(sum_prog/(i+1),var_prog/(i+1),i) ;
         pi = sum_prog/(i+1);
         Result << pi << " " << err << endl; // salvo media e sua varianza su file
      } else cerr << "PROBLEM: Unable to open results.out" << endl;
   }

   // Conto le cifre significative prima di stampare a video
   double errcut = err;
   int signfc=1;
   while(errcut < 1){
      signfc++;
      errcut *= 10;
   }
   // Stampo l'ultimo risultato
   cout << setprecision(signfc) << "Last extimation of PI: " << pi << setprecision(1) << " +/- " << err << endl;

   Result.close();
   Result2.close();

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

