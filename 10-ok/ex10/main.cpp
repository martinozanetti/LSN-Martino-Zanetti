/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex10

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "main.h"
#include "../../../../anaconda3/include/mpi.h"
//#include "/include/mpi.h" // compila ma VS me lo segna in rosso

using namespace std;

int main(int argc, char* argv[]){

// ====== INIZIALIZZAZIONE PROBLEMA ======
    Random rnd;
    rnd.SetSeed();
    rnd.SetPrimesCouple(23);
    Problemset pr;
    Population pop;
    int nPop = pop.GetNPop();    // nuomero di cromosomi (= individui nella popolazione)
    int nPas = NGeneration;      // numero di generazioni consecutive

    // genero città e popolazione di cromosomi associata, valuto e metto in ordine
    //pr.GenCircCities();       // <<< COMMENTARE ALTERNATIVAMENTE QUESTE DUE RIGHE A SECONDA DEL DESIDERIO
    //pr.GenSquareCities();     // <<< COMMENTARE ALTERNATIVAMENTE QUESTE DUE RIGHE A SECONDA DEL DESIDERIO
    pr.LoadCities("American_capitals.dat"); // <<< COMMENTARE ALTERNATIVAMENTE QUESTE DUE RIGHE A SECONDA DEL DESIDERIO
    pop.Birth(rnd);
    pr.EvalAll(pop);    
    pr.SortPop(&pop);

// ====== INIZIALIZZAZIONE MPI ======
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;

    rnd.SetPrimesCouple(rank); // inizializzo su coppie di semi diversi ciascun nodo

// ======== EVOLUZIONE ========
    bool migrate = true; // <<< COMMENTARE ALTERNATIVAMENTE QUESTE DUE RIGHE A SECONDA DEL DESIDERIO
    int tau_migr = 30; // ogni quante generazioni si ha una migrazione

    for(int k = 0; k < nPas; k++){ // ciclo sulle generazioni

        if(!migrate || (migrate && k%tau_migr!=0) || k==0){ // continenti indipendenti -> no migrazioni
            pop.Mutate(rnd);
            pr.EvalAll(pop);    
            pr.SortPop(&pop);
        }
        else if(migrate && k%tau_migr==0){ // continenti connessi -> migrazioni
            // ogni tau_migr generazioni un chr di una
            // popolazione in un *core casuale* dev'essere 
            // mandato in un'altro *core casuale*
            // --- migrare è analogo a mutare ---
            int h = 0;
            int giver = 0;
            int receiver = 1;
            int migrator[Ngenes];
            int itag=1;
            int expo = 20; //expon;

            // solo il core 0 fa le estrazioni
            if(rank == 0){
                h = (int)(nPop*(1-pow(rnd.Rannyu(),expo)))-1; // scelgo un individuo come lo avevo scelto per mutate
                giver = (int)rnd.Rannyu(0,size);             // scelgo il core di partenza
                receiver = (int)rnd.Rannyu(0,size);          // scelgo il core di arrivo... 
                while(receiver==giver){                        // ... diverso dal giver
                    receiver = (int)rnd.Rannyu(0,size);      
                }
                cout << endl << "Migration (h, giver, receiver): " << h << ", " << giver << ", " << receiver << endl;
            }

            MPI_Bcast(&h, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); // il core 0 manda a tutti gli altri
            MPI_Bcast(&giver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); 
            MPI_Bcast(&receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); 

            // se sono il giver, salvo su migrator e lo spedisco al receiver
            if(rank == giver){
                for(int i = 0; i< Ngenes; i++){
                    migrator[i] = pop.Chr[h].GetGen(i);
                }
                MPI_Send(migrator,Ngenes,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD);

            }
            // se sono il ricevitore lo salvo
            if(rank == receiver){
                MPI_Recv(migrator,Ngenes,MPI_INTEGER,giver,itag,MPI_COMM_WORLD,&stat);
                pop.Chr[h].SetGen(migrator); 
                pop.Chr[h].Check();
            }
        
            // dopo è uguale
            pr.EvalAll(pop);    
            pr.SortPop(&pop);  
        }
        
        // salvo i risultati
        pr.PrintCities(k, pop.Chr[nPop-1], rank); //scegliendo Chr[nPop-1], prendo il migliore
        pr.PrintBestsLenAve(k, nPop/2, pop, rank);
        pr.PrintBestLen(k, pop, rank);
        
        if(k%10 == 0) cout << "Generation n: " << k  << " of " << nPas << endl;
    }
    

    MPI_Finalize();
    
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