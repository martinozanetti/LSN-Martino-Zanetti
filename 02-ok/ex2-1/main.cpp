/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex2-1

#include <iostream>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include "../../mylib/mylib.h"
#include <cmath>

using namespace std;

int main (){

    int M=10000;        // number of throws in each block
    int N=100;           // Number of blocks = number of experiments
   
    double sum_prog = 0; // variabili di appoggio per la media
    double var_prog = 0;

    double sum_prog_inv = 0; // per metodo di inversione della cumulativa
    double var_prog_inv = 0;

    Random rnd;  
    rnd.SetSeed();

    ofstream WriteResult;
    WriteResult.open("results.out"); 
    ofstream WriteResultIS;
    WriteResultIS.open("results_IS.out"); 

    // Ciclo dei blocchi
    for(int i=0; i<N; i++){
        double mean = 0; 
        double mean_inv = 0; 

    // Ciclo dei lanci in un blocco: x in intervallo 0,1, poi valuto la funzione e conto+
        for(int j=0; j<M; j++){ 
            double x = rnd.Rannyu(); // metodo media
            double y = M_PI/2.*cos(M_PI*x/2.);
            mean += y;

            double x_inv = rnd.Linear(-2., 2., 0., 1.); // estraggo con densità d(x)
            y = M_PI/2.*cos(M_PI*x_inv/2.)/(2-2.*x_inv);
            mean_inv += y;
        }

    // calcolo la media dividendo per il totale dei lanci: l'area è quel risultato * la lunghezza dell'intervallo

        mean = double(mean/M);  // media
        sum_prog += mean;       // somma medie
        var_prog += mean*mean;  // somma medie quadrate

        mean_inv = double(mean_inv/M); // G_i
        sum_prog_inv += mean_inv;
        var_prog_inv += mean_inv*mean_inv;
    
        // salvo su file: media e suo errore; sigma e suo errore
        if (WriteResult.is_open()){
            WriteResult << sum_prog/(i+1)     << " " << error(sum_prog    /(i+1),var_prog    /(i+1),i) << endl; //
        } else {cerr << "PROBLEM: Unable to open results.out" << endl;}
        if (WriteResultIS.is_open()){
            WriteResultIS << sum_prog_inv/(i+1) << " " << error(sum_prog_inv/(i+1),var_prog_inv/(i+1),i) << endl; //
        } else {cerr << "PROBLEM: Unable to open results_IS.out" << endl;}
   }
   WriteResult.close();
   WriteResultIS.close();

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

