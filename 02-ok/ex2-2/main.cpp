/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex2-2

#include <iostream>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include "../../mylib/mylib.h"
#include <cmath>

using namespace std;

int main (){

    int walks_tot = 10000;              // numero di random-walks (ripetizioni)
    int steps = 100;                    // numero di passi in ogni random-walk
    int blocks = 100;                   // numero di blocchi ("raggruppamenti" di random-walks)
    int walks_per_block = (int)walks_tot/blocks;  // numero di random_walks per blocco
   
    double mean_radius_discr[steps][blocks] = {0}; // variabili di appoggio per la versione discreta
    double ave_radius_discr[steps] = {0};
    double ave2_radius_discr[steps] = {0};

    double mean_radius_cont[steps][blocks] = {0}; // variabili di appoggio per la versione continua
    double ave_radius_cont[steps] = {0};
    double ave2_radius_cont[steps] = {0};

    Random rnd;  
    rnd.SetSeed();

    ofstream WriteResultDiscr;
    WriteResultDiscr.open("discrete.out");
    ofstream WriteResultCont;
    WriteResultCont.open("continuous.out");

    // Ciclo dei blocchi: divido i miei lanci in un numero di blocchi sensato
    for(int i=0; i<blocks; i++){

    // ciclo dei random-walks (ogni volta riparto dall'origine): calcolo sqrt(mean_in_block(modulo_quadro)) in funzione dello step
        for(int j=0; j<walks_per_block; j++){

            double pos_discr[3] = {0};
            double pos_cont[3] = {0};

    // Ciclo dei passi in un walk
            for(int k=0; k<steps; k++){ 

                // RETICOLO DISCRETO
                int step_dir  = int(rnd.Rannyu(0.,3.)); // direzione dello step: 0=x, 1=y, 2=z
                int step_sign = int(rnd.Rannyu(1.,3.)); // verso dello step: 1=avanti, 2=indietro

                if      (step_sign == 1)    step_sign = +1;
                else if (step_sign == 2)    step_sign = -1;

                pos_discr[step_dir] += step_sign;

                // per ogni random-walk, salvo il modulo quadro della posizione al passo k-esimo, mediato sul numero random-walks per blocco
                mean_radius_discr[k][i] += ( pow(pos_discr[0],2) + pow(pos_discr[1],2) + pow(pos_discr[2],2) ) / walks_per_block;

                // AMBIENTE CONTINUO
                // direzioni random su una semisfera...
                double cosThetaSq = rnd.Rannyu();
                double cosTheta = sqrt(cosThetaSq);
                double senTheta = sqrt(1.0 - cosThetaSq);
                double phi   = rnd.Rannyu(0., 2.*M_PI);
                // estensione al secondo emisfero
                int up_down = ((int)rnd.Rannyu(0.,2.)*2-1); // +-1

                pos_cont[0] += senTheta*cos(phi);
                pos_cont[1] += senTheta*sin(phi);
                pos_cont[2] += cosTheta*up_down;

                /* // codice precedente: non restituisce una distribuzione uniforme
                double theta = rnd.Rannyu(0., M_PI);
                double phi   = rnd.Rannyu(0., 2.*M_PI);

                pos_cont[0] += sin(theta)*cos(phi);
                pos_cont[1] += sin(theta)*sin(phi);
                pos_cont[2] += cos(theta);
                */

                // per ogni random-walk, salvo il modulo quadro della posizione al passo k-esimo, mediato sul numero random-walks per blocco
                mean_radius_cont[k][i] += ( pow(pos_cont[0],2) + pow(pos_cont[1],2) + pow(pos_cont[2],2) ) / walks_per_block;

            }    
        }
    }

    // ho finito di riempire la matrice delle radici delle medie nei blocchi delle distanze quadrate

    if (WriteResultDiscr.is_open() && WriteResultCont.is_open()){

        for (int i = 0; i < steps; i++)
        {
            // fissato lo step i-esimo calcolo media e media quadratica sui blocchi
            for (int j = 0; j < blocks; j++)
            {
                double appo = mean_radius_discr[i][j];
                ave_radius_discr[i] += appo/blocks;
                ave2_radius_discr[i] += appo*appo/blocks;

                double appo_cont = mean_radius_cont[i][j];
                //cout << appo_cont << endl;
                ave_radius_cont[i] += appo_cont/blocks;
                ave2_radius_cont[i] += appo_cont*appo_cont/blocks;
            
            }
            
            double err = 0;
            if(i == 0) WriteResultDiscr << sqrt(ave_radius_discr[i]) << " " << 0 << endl;
            else
            {
                // propagazione errore: normalizzo per 1/2\sqrt{ave}
                err = error(ave_radius_discr[i], ave2_radius_discr[i], blocks)/(2*sqrt(ave_radius_discr[i]));
                WriteResultDiscr << sqrt(ave_radius_discr[i]) << " " << err << endl;
            }
            
            if(i == 0) WriteResultCont << sqrt(ave_radius_cont[i]) << " " << 0 << endl;
            else
            {
                // propagazione errore: normalizzo per 1/2\sqrt{ave}
                err = error(ave_radius_cont[i], ave2_radius_cont[i], blocks)/(2*sqrt(ave_radius_cont[i]));
                WriteResultCont << sqrt(ave_radius_cont[i]) << " " << err << endl;
            }

        }

    } else {
        if (!WriteResultDiscr.is_open()) cerr << "PROBLEM: Unable to open discrete.out" << endl;
        if (!WriteResultCont.is_open()) cerr << "PROBLEM: Unable to open continuous.out" << endl;
    }
   
   WriteResultDiscr.close();
   WriteResultCont.close();

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

