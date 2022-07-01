/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex3-1

#include <iostream>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include "../../mylib/mylib.h"
#include <cmath>
#include <iomanip>


using namespace std;

int main (){

    int M = 10000;                    // numero totale di ripetizioni
    int blocks = 100;                   // numero di blocchi
    int steps = (int)M/blocks;          // numero di step in ogni blocco

    Random rnd;  
    rnd.SetSeed();

    ofstream WriteResultDirect;
    WriteResultDirect.open("direct.out");
    ofstream WriteResultDiscret;
    WriteResultDiscret.open("discretized.out");


    double S = 0; // variabile di appoggio
    double S_0 = 100; // asset price al tempo 0
    double T = 1; // expiry time
    double k = 100; // prezzo concordato
    double r = 0.1; // risk-free interest rate
    double sigma = 0.25; // volatilità

    // variabili per diretto
    double C = 0; // call-option price
    double var_prog_C = 0; // errori statistici
    double mean_prog_C = 0;

    double P = 0; // put-option price
    double var_prog_P = 0; // errori statistici
    double mean_prog_P = 0;

    // variabili per discreto
    double C_d = 0; // call-option price
    double var_prog_C_d = 0; // errori statistici
    double mean_prog_C_d = 0;

    double P_d = 0; // put-option price
    double var_prog_P_d = 0; // errori statistici
    double mean_prog_P_d = 0;

    // ciclo diretto: genero direttamente la media al tempo T del costo della mia opzione...
    for(int i = 0; i < blocks; i++)
    {
        C = 0;
        P = 0;
        C_d = 0;
        P_d = 0;

        // ciclo nel singolo blocco
        for (int j = 0; j < steps; j++)
        {
            double z = rnd.Gauss(0,1);                              
            S = S_0*exp((r-sigma*sigma/2.)*T + sigma*z*pow(T,0.5)); // gennero il valore qui se non ho i blocchi...
            C += exp(-r*T)*max(0.,S-k)/steps;                       // esponenzio il valore, cioè lo sconto al presente
            P += exp(-r*T)*-min(0.,S-k)/steps;
            //P += exp(-r*T)*max(0.,k-S)/steps;


            //... altrimenti entro nel ciclo della discretizzazione
            int time_steps = 100;
            for (int h = 0; h < time_steps; h++)
            {
                double dt = T/100.; // passettino discreto di tempo
                z = rnd.Gauss(0,1);
                if (h == 0) S = S_0*exp((r-sigma*sigma/2.)*dt + sigma*z*pow(dt,0.5)); // per il primo passo ho S_1(S_0)...
                else        S = S  *exp((r-sigma*sigma/2.)*dt + sigma*z*pow(dt,0.5)); // ... per gli altri ho S_h+1(S_h)
            }
            C_d += exp(-r*T)*max(0.,S-k)/steps;
            P_d += exp(-r*T)*-min(0.,S-k)/steps;
            //P += exp(-r*T)*max(0.,k-S)/steps;

        } 
        
        // metodo diretto
        mean_prog_C += C;
        var_prog_C += C*C;

        mean_prog_P += P;
        var_prog_P += P*P;

        // metodo discretizzato
        mean_prog_C_d += C_d;
        var_prog_C_d += C_d*C_d;

        mean_prog_P_d += P_d;
        var_prog_P_d += P_d*P_d;

        // salvo su file
        if (WriteResultDirect.is_open() && WriteResultDiscret.is_open())
        {            
            if(i == 0) WriteResultDirect << mean_prog_C/(i+1) << " " << 0 << " " << mean_prog_P/(i+1) << " " << 0 << endl;
            else WriteResultDirect << mean_prog_C/(i+1) << " " << error(mean_prog_C/(i+1), var_prog_C/(i+1), i) << " " //
                                   << mean_prog_P/(i+1) << " " << error(mean_prog_P/(i+1), var_prog_P/(i+1), i) << " " << endl;
            
            if(i == 0) WriteResultDiscret << mean_prog_C_d/(i+1) << " " << 0 << " " << mean_prog_P_d/(i+1) << " " << 0 << endl;
            else WriteResultDiscret << mean_prog_C_d/(i+1) << " " << error(mean_prog_C_d/(i+1), var_prog_C_d/(i+1), i) << " " //
                                    << mean_prog_P_d/(i+1) << " " << error(mean_prog_P_d/(i+1), var_prog_P_d/(i+1), i) << " " << endl;

        } else {
            if (!WriteResultDirect.is_open()) cerr << "PROBLEM: Unable to open direct.out" << endl;
            if (!WriteResultDiscret.is_open()) cerr << "PROBLEM: Unable to open discretized.out" << endl;
        }
    }

    // medie finali = tra loro e al risultato analitico?
    mean_prog_C = mean_prog_C/blocks; 
    mean_prog_P = mean_prog_P/blocks; 

    mean_prog_C_d = mean_prog_C_d/blocks; 
    mean_prog_P_d = mean_prog_P_d/blocks; 

    cout << endl << "Confronto risultati numerici con risultati analitici (direct, discretized): " << endl; 
    cout << "Call (an: 14.98): [" << setprecision(3) << mean_prog_C << " +/- "   << setprecision(2) << error(mean_prog_C/(blocks+1),   var_prog_C/(blocks+1),   blocks) << "];  [" //
                                  << setprecision(3) << mean_prog_C_d << " +/- " << setprecision(2) << error(mean_prog_C_d/(blocks+1), var_prog_C_d/(blocks+1), blocks) << "]" << endl;
    cout << "Put  (an: 5.46):  [" << setprecision(3) << mean_prog_P << " +/- "   << setprecision(2) << error(mean_prog_P/(blocks+1),   var_prog_P/(blocks+1),   blocks) << "]; ["  //
                                  << setprecision(3) << mean_prog_P_d << " +/- " << setprecision(2) << error(mean_prog_P_d/(blocks+1), var_prog_P_d/(blocks+1), blocks) << "]" << endl;
   
    WriteResultDirect.close();
    WriteResultDiscret.close();

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

