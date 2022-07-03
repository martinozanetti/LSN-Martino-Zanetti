
#include <iostream>
#include <fstream>
#include "../random/random.h"
#include "parameters.h"

//========================================================//
//       CHROMOSOME                                       //
//========================================================//

#ifndef __CHROMOSOME__
#define __CHROMOSOME__

class Chromosome {

private:
  static const int Ng = Ngenes; // regolare da parameters
  int Gen[Ng] = {0};

public:
  double Ftn = 0;
  double Len = 0;
  const double *Pm = pmut; // basic mutation probability

  // constructors
  Chromosome();

  // destructor
  ~Chromosome();

  Chromosome& operator= (const Chromosome& chr);

  // methods

  // set e get
  void SetGen(int vec[Ng]);
  void SetFitness(double f);
  void SetLen(double l);
  double GetFtn();
  double GetLen();
  int  GetNg();
  int  GetGen(int i);

  // altri
  void Fill(Random &rnd);
  void Empty();
  void Print();

  // mutations
  void PairPermut(int n, int m);
  void Shift(int pos, int m, int n);
  void Shift2(int shift);
  void MPermut(int pos1, int pos2, int m);
  void Inversion(int pos, int m);
  void Crossover(int pos, int len, Chromosome parent2);
  // accessori alle mutazioni
  int Pbc(int pos);
  bool Check();


};

#endif // __CHROMOSOME__

//========================================================//
//       POPULATION                                       //
//========================================================//

#ifndef __POPULATION__
#define __POPULATION__

class Population {

private:
  static const int Npop = NindPop; // regolare da parameters

public:
  double BestLen;
  double BestsLenAve;
  Chromosome Chr[Npop];

  // constructors
  Population();

  // destructor
  ~Population();

  // methods
  int GetNPop();
  void Birth(Random &rnd);
  void Mutate(Random &rnd);
};

#endif // __POPULATION__


//========================================================//
//       PROBLEMSET                                       //
//========================================================//

#ifndef __PROBLEMSET__
#define __PROBLEMSET__

class Problemset {

  private:
    static const int Ncit = Ngenes; // regolare da parameters
    double Xcit[Ncit];
    double Ycit[Ncit];

  public:
    Problemset(); // costruttore
    ~Problemset(); // distruttore

    void GenCircCities();
    void GenSquareCities();
    void PrintCities(int generation, Chromosome chr, int rank);
    void PrintBestsLenAve(int generation, int part, Population pop, int rank);
    void PrintBestLen(int generation, Population pop, int rank);
    void EvalFitness(Chromosome &chr);
    void EvalAll(Population &pop);
    void SortPop(Population *pop);

};

#endif // __PROBLEMSET__

//========================================================//
//       OTHER FUNCTIONS                                  //
//========================================================//

int  partition(Population *pop, int low, int high);
void quickSort(Population *pop, int low, int high);