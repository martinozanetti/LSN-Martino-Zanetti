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
    pr.GenSquareCities();     // <<< COMMENTARE ALTERNATIVAMENTE QUESTE DUE RIGHE A SECONDA DEL DESIDERIO
    pop.Birth(rnd);
    pr.EvalAll(pop);    
    pr.SortPop(&pop);

// ====== INIZIALIZZAZIONE MPI ======
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    rnd.SetPrimesCouple(rank); // inizializzo su coppie di semi diversi ciascun nodo

// ======== EVOLUZIONE ========
    // ciclo sulle generazioni
    for(int k = 0; k < nPas; k++){
        
        pop.Mutate(rnd);
        pr.EvalAll(pop);    
        pr.SortPop(&pop);

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