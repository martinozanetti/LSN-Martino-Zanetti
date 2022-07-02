/****************************************************************
*****************************************************************
    _/    _/ _/_/_/_/      Numerical Simulation Laboratory
   _/_/_/_/      _/       Physics Department
  _/ _/ _/    _/         Universita' degli Studi di Milano
 _/    _/  _/           Martino Zanetti
_/    _/ _/_/_/_/      email: martino.zanetti@gmail.com
*****************************************************************
*****************************************************************/
// ex9

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "main.h"

using namespace std;

int main (){

// ====== INIZIALIZZAZIONE ======
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

// ======== EVOLUZIONE ========
    // ciclo sulle generazioni
    for(int k = 0; k < nPas; k++){
        
        pop.Mutate(rnd);
        pr.EvalAll(pop);    
        pr.SortPop(&pop);

        // salvo i risultati
        pr.PrintCities(k, pop.Chr[nPop-1]); //scegliendo Chr[nPop-1], prendo il migliore
        pr.PrintBestsLenAve(k, nPop/2, pop);
        pr.PrintBestLen(k, pop);

        if(k%10 == 0) cout << "Generation n: " << k  << " of " << nPas << endl;
    }

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