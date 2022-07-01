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

#include <iostream>
#include <fstream>
#include <cmath>
//#include <string>
#include <iomanip>
#include "main.h"

//using namespace std;

//==============================================

int main (){

    Input();

    //equilibrate
    for(int i = 0; i<1000; i++)
    {
        Move();
    }

    for(int i = 0; i < nblk; i++) // ciclo sui blocchi
    {     
        ResetAll();
        
        for (int j = 0; j < nstep; j++) // ciclo nel blocco
        {           
            SavePos();  // current position     
            Move();
            Accumulate();
        }            

        PrintAccRate(i); // acceptance rate
        BlockAverages();
        SaveDist(i); // save out block-averaged distances

    }

    return 0;
}

// ============== Functions ===================

void Input(void)
{
  ifstream ReadInput;

  cout << ",-------------------------------------," << endl;
  cout << "| Hydrogen atom                       |" << endl;
  cout << "| Monte Carlo simulation (Metropolis) |" << endl;
  cout << "'-------------------------------------'" << endl;

//Read seed for random numbers
   rnd.SetSeed();
  
//Read input informations
  ReadInput.open("input.dat");
    
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> L;
  ReadInput >> r[0];
  ReadInput >> re[0];
  ReadInput >> gauss;
  //ReadInput >> np_rend;

  if (!gauss) tr_pr = "unif";
  else if (gauss) tr_pr = "gauss";

  steps = nstep*nblk;
  
  //mod = (int)steps/np_rend;
  mod=1;

  cout << "Total number of steps = " << steps << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  cout << "Step lenght = " << L << endl;
  cout << "Initial position (GS) = " << r[0] << endl;
  cout << "Initial positin (ExSt) = " << re[0] << endl;
  //cout << "Num punti da renderizzare = " << np_rend << endl;
  cout << endl;

  ReadInput.close();

}

void ResetAll(){

    d = 0; // distanze dal centro
    de = 0; 
    attempted_gs = 0; // contatori acceptance rate
    attempted_es = 0;
    accepted_gs = 0;
    accepted_es = 0;
}

void SavePos(){

    ofstream WritePosGS;
    WritePosGS.open("GS/pos_"+ tr_pr +".out", ios::app);
    ofstream WritePosES;
    WritePosES.open("ES/pos_"+ tr_pr +".out", ios::app);

    // salvo la posizione attuale
    for(int k = 0; k<3; k++){
        y[k] = r[k]; // GS
        ye[k] = re[k]; // ES
    
        // anche su file ogni mod passi
        if(contastep%mod == 0) 
        {
            WritePosGS << y[k] << " ";
            WritePosES << ye[k] << " "; 
        }
    }
    if(contastep%mod == 0)
    {
        WritePosGS << endl;
        WritePosES << endl;
    }

    WritePosGS.close();
    WritePosES.close();
}

void ProposeStep(){

    /*
    //                    ^
    //                    |
    //            ________|________
    //           |        |        |
    //-----------,========|========,------------->
    //          -L                +L
    //
    */
    if(!gauss){
        for(int k = 0; k<3; k++){
            dr_g[k] = {rnd.Rannyu(-L,L)}; // uniform transition probability in un cubo di lato 2L
            dr_e[k] = {rnd.Rannyu(-L,L)}; // idem per stato eccitato
        }
    }
    /*
    //                    ^
    //                    |
    //                    _
    //                  / | \
    //               /    |    \
    //             /______|______\      2L=FWHM
    //           /        |        \
    //        /  |        |        |  \
    //-----<-----,========|========,----->-------->
    //         -sigma            +sigma
    //                
    //                 sigma = L
    //
    */
    else if(gauss){
        for(int k = 0; k<3; k++){
            dr_g[k] = {rnd.Gauss(0,L)}; // gaussian transition probability con sigma = L
            dr_e[k] = {rnd.Gauss(0,L)}; // idem per stato eccitato
        }
    }
}

void Move(){

    ProposeStep();

    // provo a spostarmi
    contastep++;
    for(int k = 0; k<3; k++){
        x[k] = y[k] + dr_g[k]; // GS
        xe[k] = ye[k] + dr_e[k]; // ES
    }

    // calcolo probabilità e accetto secondo metropolis

    //GS
    attempted_gs++;
    q = prob_gs(x)/prob_gs(y);
    A = min(1,q);

    if (A==1) // accetto direttamente
    {
        for(int k = 0; k<3; k++) r[k] = x[k];
        accepted_gs++;
    }
    else // accetto con probabilità A
    {
        if(rnd.Rannyu() < A)
        {
            for(int k = 0; k<3; k++) r[k] = x[k];
            accepted_gs++;
        }
    }

    // ES
    attempted_es++;
    qe = prob_exc(xe)/prob_exc(ye);
    Ae = min(1,qe);

    if (Ae==1) // accetto direttamente
    {
        for(int k = 0; k<3; k++) re[k] = xe[k];
        accepted_es++;
    }
    else // accetto con probabilità A
    {
        if(rnd.Rannyu() < Ae)
        {
            for(int k = 0; k<3; k++) re[k] = xe[k];
            accepted_es++;
        }
    }
}

void Accumulate(){
    // accumulo nel blocco
    d += sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    de += sqrt(re[0]*re[0]+re[1]*re[1]+re[2]*re[2]);
}

void PrintAccRate(int i){
    // stampo acceptance rate
    cout << "Block # " << i+1 << endl;
    cout << "                    (GrSt)  (ExSt)  " << endl;
    cout << "Acceptance rates:   " << (double)accepted_gs/attempted_gs << "  " << (double)accepted_es/attempted_es << endl;
    cout << "-----------------------------------" << endl;
}

void BlockAverages(){
    // medie di blocco
    d /= nstep;   
    mean_prog_d += d;
    var_prog_d += d*d;

    de /= nstep; 
    mean_prog_de += de;
    var_prog_de += de*de;
}

void SaveDist(int i){

    ofstream WriteResultGS;
    WriteResultGS.open("GS/dist_"+ tr_pr +".out", ios::app);
    ofstream WriteResultES;
    WriteResultES.open("ES/dist_"+ tr_pr +".out", ios::app);

    // salvo su file
    // GS
    if (WriteResultGS.is_open())
    {            
        if(i == 0) WriteResultGS << d << " " << mean_prog_d/(i+1) << " " << 0 << endl;
        else WriteResultGS << d << " " << mean_prog_d/(i+1) << " " << error(mean_prog_d/(i+1), var_prog_d/(i+1), i) << " " << endl;
        
    } else {
        if (!WriteResultGS.is_open()) cerr << "PROBLEM: Unable to open resultGS.out" << endl;
    }
    // ES
    if (WriteResultES.is_open())
    {            
        if(i == 0) WriteResultES << de << " " << mean_prog_de/(i+1) << " " << 0 << endl;
        else WriteResultES << de << " " << mean_prog_de/(i+1) << " " << error(mean_prog_de/(i+1), var_prog_de/(i+1), i) << " " << endl;
        
    } else {
        if (!WriteResultES.is_open()) cerr << "PROBLEM: Unable to open resultES.out" << endl;
    }    

    WriteResultGS.close();
    WriteResultES.close();
}

double prob_gs(double x[3])
{
    double d = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    double psi = pow(M_E,-d)/sqrt(M_PI);

    return psi*psi;
}

double prob_exc(double x[3])
{
    double d = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    double costheta = x[2]/d;

    double psi = 1./8.*sqrt(2./M_PI)*d*pow(M_E,-d/2)*costheta;

    return psi*psi;
}

double min(double s, double t)
{
    if (s <= t) return s;
    else return t;
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