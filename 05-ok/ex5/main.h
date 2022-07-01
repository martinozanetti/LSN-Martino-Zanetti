/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex5

#ifndef __EX5__
#define __EX5__

//Random numbers & my library
#include "../../random/random.h"
#include "../../mylib/mylib.h"
#include <string>

using namespace std;

Random rnd;  

int nstep = 1e4;                      // numero di step per blocco
int nblk = 100;                   // numero di blocchi
int steps = nstep*nblk;          // numero di step totale
double L = 1.5;                   // lato del cubo visitabile / sigma per caso gaussiano
bool gauss;
int np_rend = 10000;
//int mod = (int)steps/np_rend;
int mod = 1;
int contastep = 0;

string tr_pr;

// variabili di appoggio per Ground State
double d = 0;           // distanza dal centro
double var_prog_d = 0;  // appoggio per errori statistici
double mean_prog_d = 0;

double r[3] = {10.0,0.0,0.0}; // posizione iniziale
double y[3];                 // posizione prima del passo
double x[3];                 // posizione dopo il passo

double dr_g[3];

double q, A;
int attempted_gs = 0, accepted_gs = 0;

// variabili di appoggio per Excited State
double de = 0;          // distanza dal centro
double var_prog_de = 0; // appoggio per errori statistici
double mean_prog_de = 0;

double re[3] = {20.0,0.0,0.0}; // posizione iniziale
double ye[3];                 // posizione prima del passo
double xe[3];                 // posizione dopo il passo

double dr_e[3];

double qe, Ae;
int attempted_es = 0, accepted_es = 0;

//functions
void Input(void);
double prob_gs(double*);
double prob_exc(double*);
double min(double, double);
void ResetAll();
void SavePos();
void Move();
void Accumulate();
void PrintAccRate(int blknum);
void BlockAverages();
void SaveDist(int blknum);


#endif

/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/