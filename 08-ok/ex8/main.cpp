/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex8

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "main.h"

using namespace std;

ofstream WriteResult;
ofstream WritePos;

ofstream WriteTraj;

//==============================================

int main (){

    Input();

    if(!SA)
    {
        MeasureAverageH();
    }

    if(SA)
    {
        WriteTraj.open("traj.out");

        /* ============================================================
        |
        | Partire con mu, sigma, T(*) sensati. (Dobbiamo scegliere 
        | deltaT, dim passo, ...)
        | Equilibrare.
        | Poi fare x passi (che può essere anche uno solo!), 
        | abbassare T, fare altri passi, ecc... 
        | ...fino a temperature bassissime dove ci fermiamo.
        | 
        |
        | Fare un ulteriore metropolis per prendere o no 
        | i nuovi parametri: 
        | se l'energia(**) diminuisce allora deltaH ci fa accettare, 
        | altrimenti accettiamo con probabilità q.
        |
        |  (*) Nota che T è un parametro fittizio, che deve misurare la 
        |      capacità di andare a visitare stati (mu, sgm) lontani 
        |      dall'optimum.
        |      Riflettere su come implementarlo! 
        |    
        | (**) val aspett di H, cioè la H media, con rispettivo errore. 
        |      Nota che qui il deltaH deve essere maggiore di errH, 
        |      altrimenti non ha senso procedere.
        |      In effetti, quando diventano confrontabili, 
        |      è meglio fermarsi.
        |
        |============================================================== */

        // prendo mu, sigma, T (>> beta) iniziali sensate
            // (si impostano da input.dat)

        // equilibrare facendo qualche passo (soprattutto se sono lontano da un punto di minimo)
        // se parto da dentro la buca, una manciata (5) di passi dovrebbero essere suff.

        int eq = 10;
        for(int i = 0; i < eq; i++)
        {
            MetroMove();
            cout << "equilibrio: " << i+1 << "/" << eq << endl << endl; 
        }
        cout << "==================" << endl << endl;

        H = 0; // riazzero tutto prima di iniziare SA
        attempted = 0;
        accepted = 0;

        // ciclo sul beta crescente (raffreddamento) con un certo numero di passi per ogni beta, contando gli step totali.
        // nel ciclo ogni volta calcolo H con rispettivo errore.
        // ogni passo deve essere giudicato alla metropolis: se H diminuisce accetto, altrimenti accetto con prob q.


        while (2*Lmu/beta > errMu or 2*Lsgm/beta > errSgm)    // mi fermo quando a una data temperatura le oscillazioni percentuali di H sono minori di un tot
        {            
            for(int i = 0; i < step_in_beta; i++)
            {
                oldH = mean_prog_H/nblk;

                // equivalente di passo tentativo: esploro mu e sgm in funzione di beta
                d_mu = rnd.Rannyu(-Lmu,Lmu)/beta; // a basse T >> alte beta, l'esplorazione è meno ampia.
                d_sgm = rnd.Rannyu(-Lsgm,Lsgm)/beta;

                mu += d_mu;
                sgm += d_sgm;

                MeasureAverageH();
                attempted++;
                
                // calcolo probabilità
                newH = mean_prog_H/nblk;
                errH = error(mean_prog_H/nblk,
                             var_prog_H/nblk, 
                             (nblk-1));
                deltaH = oldH - newH;

                q = exp(beta*(deltaH)); // con +, perchè devo andare verso il minimo di H
                A = min(1,q);

                // giudico con metropolis, confrontantdo newH con oldH
                if (A==1) // accetto direttamente
                {
                    oldH = newH;
                    accepted++;
                }
                else // accetto con probabilità A
                {
                    if(rnd.Rannyu() < A)
                    {
                        oldH = newH;
                        accepted++;
                    } else {mu -= d_mu; sgm -= d_sgm;} // se non va, ripristino i parametri
                }
            
                step++;

            }

            cout << "step:    " << step << endl 
                 << "beta:    " << beta << endl 
                 << "H:       " << oldH << endl
                 << "errH:    " << errH << endl
                 //<< "errH(%): " << errH/abs(oldH)*100 << endl
                 << "mu:      " << mu << endl
                 << "errMu:   " << 2*Lmu/beta << endl
                 << "sgm:     " << sgm << endl
                 << "errSgm:  " << 2*Lsgm/beta << endl
                 << endl; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ATTENZIONE QUESTO RISULTATO NON È COERENTE CON LA MIA RICHIESTA PERCENTUALE
            WriteTraj << beta << " " << mu << " " << sgm << " " << oldH  << " " << errH << endl;

            // abbasso temperatura
            beta += db;
            
            // all'ultima temperatura raggiunta, mi fermo e mi muovo un po in giro per capire il vero errore di mu e sigma
            if (2*Lmu/beta <= errMu and 2*Lsgm/beta <= errSgm){

                cout << "Mi muovo un po' fissata la temperatura più bassa, per valutare la STD di mu e sigma." << endl;
                
                ofstream WriteTraj_last;
                WriteTraj_last.open("traj_last.out");


                for(int i = 0; i < step_in_beta*amp; i++)
                {
                    if (i%1000 == 0) cout << i << "/" << step_in_beta*amp << endl;

                    oldH = mean_prog_H/nblk;

                    // equivalente di passo tentativo: esploro mu e sgm in funzione di beta
                    d_mu = rnd.Rannyu(-Lmu,Lmu)/beta; // a basse T >> alte beta, l'esplorazione è meno ampia.
                    d_sgm = rnd.Rannyu(-Lsgm,Lsgm)/beta;

                    mu += d_mu;
                    sgm += d_sgm;

                    MeasureAverageH();
                    attempted++;
                
                    // calcolo probabilità
                    newH = mean_prog_H/nblk;
                    errH = error(mean_prog_H/nblk,
                             var_prog_H/nblk, 
                             (nblk-1));
                    deltaH = oldH - newH;

                    q = exp(beta*(deltaH)); // con +, perchè devo andare verso il minimo di H
                    A = min(1,q);

                    // giudico con metropolis, confrontantdo newH con oldH
                    if (A==1) // accetto direttamente
                    {
                        oldH = newH;
                        accepted++;
                    }
                    else // accetto con probabilità A
                    {
                        if(rnd.Rannyu() < A)
                        {
                            oldH = newH;
                            accepted++;
                        } else {mu -= d_mu; sgm -= d_sgm;} // se non va, ripristino i parametri
                    }
                    step++;

                    WriteTraj_last << i << " " << mu << " " << sgm << " " << oldH << endl;

                    // calcolo errore su sgm: <<< attenzione ragionare sulla correlazione
                    mean_sgm += sgm;
                    var_sgm += sgm*sgm;
                    mean_mu += mu;
                    var_mu += mu*mu;
                }

                WriteTraj_last.close();

                mean_sgm /= step_in_beta*amp;
                var_sgm /= step_in_beta*amp;
                mean_mu /= step_in_beta*amp;
                var_mu /= step_in_beta*amp;
            }

            
        }
        cout << endl
            << ",=======================================" << endl
            << "| Minimo di H: " << oldH << " +/- " << errH << endl
            << "| mu:          " << mean_mu <<   " +/- " << error(mean_mu, var_mu, step_in_beta*amp) << endl
            << "| sigma:       " << mean_sgm  << " +/- " << error(mean_sgm,var_sgm,step_in_beta*amp) << endl
            << "'=======================================" << endl << endl;
    }
    WriteTraj.close();

    return 0;
}

// ================================================================================

// Functions
void MeasureAverageH()
{
    WriteResult.open("result.out");
    WritePos.open("pos.out");

    mean_prog_H = 0;
    var_prog_H = 0;

    // ciclo sui blocchi
    for(int i = 0; i < nblk; i++)
    {     
        H = 0; // azzero H prima di iniziare il blocco
        attempted = 0;
        accepted = 0;

        // cliclo nel singolo blocco
        for (int j = 0; j < steps; j++)
        {
            MetroMove();
            
            // accumulo H nel blocco
            H += EvalH(y);
            
            // salvo dove sono
            if(!SA && j%(steps/histofill_blk) == 0) WritePos << y << endl;
            
        }

        // stampo acceptance rate
        if(!SA)
        {
            cout << "Block # " << i+1 << endl;
            cout << "Acceptance rate:   " << (double)accepted/attempted << endl;
            cout << "-----------------------------------" << endl;
        }
                    

        // medie di blocco
        H /= steps;   
        mean_prog_H += H;
        var_prog_H += H*H;

        // salvo su file
        if(!SA)
        {
        if (WriteResult.is_open())
        {            
            if(i == 0) WriteResult << H << " " << mean_prog_H/(i+1) << " " << 0 << endl;
            else WriteResult << H << " " << mean_prog_H/(i+1) << " " << error(mean_prog_H/(i+1), var_prog_H/(i+1), i) << " " << endl;
        } else {
            if (!WriteResult.is_open()) cerr << "PROBLEM: Unable to open result.out" << endl;
        }  
        }
    }

    WriteResult.close();
    WritePos.close();
}

/*
void FileStamp(ofstream Write, int k, double val, double mean_prog_val, double var_prog_val)
{
    if (Write.is_open())
    {            
        if(k == 0) Write << val << " " << mean_prog_val/(k+1) << " " << 0 << endl;
        else Write << val << " " << mean_prog_val/(k+1) << " " << error(mean_prog_val/(k+1), var_prog_val/(k+1), k) << " " << endl;
    } else {
        if (!Write.is_open()) cerr << "PROBLEM: Unable to open output file" << endl;
    }  
}
*/

void MetroMove()
{
    dx = {rnd.Rannyu(-span,span)}; // uniform transition probability in un cubo di lato L
            
    // provo a spostarmi (x = posizione nuova)
    x = y + dx;
    attempted++;

    // calcolo prob e accetto secondo metropolis
    q = Psi2(x)/Psi2(y);
    A = min(1,q);

    if (A==1) // accetto direttamente
    {
        y = x;
        accepted++;
    }
    else // accetto con probabilità A
    {
        if(rnd.Rannyu() < A)
        {
            y = x;
            accepted++;
        }
    }
}

double EvalH(double x)
{
    double exp1 = exp(-(x-mu)*(x-mu)/(2*sgm*sgm)); 
    double exp2 = exp(-(x+mu)*(x+mu)/(2*sgm*sgm));

    double psi = exp1+exp2;
    double V = pow(x,4) - 2.5*pow(x,2);

    //double Hpsi_pot = V*psi;
    double Hpsi_kin = -0.5*( exp1*(pow(x-mu,2)/pow(sgm,4)-(pow(sgm,-2))) 
                           + exp2*(pow(x+mu,2)/pow(sgm,4)-(pow(sgm,-2)))
                           );

    return Hpsi_kin/psi + V; // << così è un filo più veloce, anche se meno chiaro

    //double Hpsi = Hpsi_kin + Hpsi_pot; 
    //return Hpsi/psi;
}

double Psi2(double x)
{
    //double exp1 = exp(-(x-mu)*(x-mu)/(2*sgm*sgm)); 
    //double exp2 = exp(-(x+mu)*(x+mu)/(2*sgm*sgm));

    //return (exp1+exp2)*(exp1+exp2);
    return pow(exp(-(x-mu)*(x-mu)/(2*sgm*sgm))+exp(-(x+mu)*(x+mu)/(2*sgm*sgm)),2);
    
}

void Input(void)
{
  ifstream ReadInput;

  cout << endl
       << ",=====================================," << endl
       << "| Simulated Annealing                 |" << endl
       << "| Monte Carlo simulation (Metropolis) |" << endl
       << "'====================================='" << endl << endl;
    
//Read input informations
  ReadInput.open("input.dat");
    
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> L;
  ReadInput >> x0;
  ReadInput >> mu;
  ReadInput >> sgm;
  ReadInput >> SA;
  ReadInput >> histofill_blk;

  ReadInput >> beta;       
  ReadInput >> db;      
  ReadInput >> step_in_beta; 

  ReadInput >> Lmu;     // larghezza iniziale passi mu
  ReadInput >> Lsgm;    // idem per sigma
  ReadInput >> errMu; // errore con cui voglio determinare mu
  ReadInput >> errSgm; // errore con cui voglio determinare mu

  ReadInput >> manualseed;
  int pr;
  ReadInput >> pr;

  if(!manualseed) rnd.SetSeed();
  else if(manualseed) {rnd.SetSeed(); rnd.SetPrimesCouple(pr);}
  
  x = x0;
  steps = nstep/nblk;
  span = L/2.0;

  cout << "Total number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << steps << endl << endl;
  cout << "Step lenght = " << L << endl;
  cout << "Initial position = " << x0 << endl;
  cout << "Num punti per riempire histo = " << histofill_blk*nblk << endl;

  ReadInput.close();

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