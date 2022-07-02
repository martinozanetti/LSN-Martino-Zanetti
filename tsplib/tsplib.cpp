
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "tsplib.h"


using namespace std;

//========================================================//
//       POPULATION                                       //
//========================================================//

/// Costruttore
Population :: Population(){}
/// Distruttore
Population :: ~Population(){}

// methods

/// riempie il membro array di cromosomi in modo randomico
void Population :: Birth(Random &rnd){
   Chromosome chr;
   for(int i =0; i < Npop; i++){
      chr.Empty();
      chr.Fill(rnd);
      Chr[i] = chr;
      Chr[i].Check();
   }
}

int Population :: GetNPop(){
   return Npop;
}

/// Metodo complesso: travasa la popolazione in una generazione successiva,
/// generata per mutazioni dei migliori individui della generazione attuale.
void Population :: Mutate(Random &rnd){

   // variabili
   Chromosome chr_appo[Npop];
   int ng = chr_appo[0].GetNg();
   const double *mp = chr_appo[0].Pm;
   int a,b,c;
   double expo = expon; // <<<<< questo parametro è scelto arbitrariamente
   int h = 0; 
   double len_av = 0;
   
   // metto la parte migliore della generazione attuale nel cromosoma di appoggio
   for (int i = 0; i < Npop; i++) {
      
      // (prima salvo lunghezze migliori)
      if(i>Npop/2) len_av += this->Chr[i].Len;          // media migliori
      if(i == Npop-1) this->BestLen = this->Chr[i].Len; // il migliore

      // estraggo uno dei migliori secondo legge di potenza, poi lo salvo
      h = (int)(Npop*(1-pow(rnd.Rannyu(),expo)))-1; // estrae un indice alto
      if(h>=Npop) h = Npop-1;                       // lim_sup
      if(h<0) h = 0;                                // lim_inf
      chr_appo[i] = this->Chr[h];
      //chr_appo[i].Check();
   }
   this->BestsLenAve = len_av/(double)(Npop/2); 

   // lo faccio evolvere e salvo per la successiva generazione
   for(int i = 0; i < Npop; i++){

      double prob = rnd.Rannyu();
      if(prob < mp[0]) {
         //cout << "--PP--" << endl;
         a = rnd.Rannyu(1,ng);
         b = rnd.Rannyu(1,ng);
         chr_appo[i].PairPermut(a,b);
      }

      prob = rnd.Rannyu();
      if(prob < mp[1]){
         //cout << "--Inv--" << endl;
         a = rnd.Rannyu(1,ng);
         b = rnd.Rannyu(1,ng);
         chr_appo[i].Inversion(a,b);
      }

      prob = rnd.Rannyu();
      if(prob < mp[2]){
         //cout << "--Shift--" << endl;
         a = rnd.Rannyu(1,ng);
         b = rnd.Rannyu(1,ng-1);
         c = rnd.Rannyu(1,ng);
         chr_appo[i].Shift(a,b,c);
      }

      prob = rnd.Rannyu();
      if(prob < mp[3]){
         //cout << "--Shift2--" << endl;
         a = rnd.Rannyu(1,ng);
         chr_appo[i].Shift2(a); 
      }

      prob = rnd.Rannyu();
      if(prob < mp[4]){
         //cout << "--MP--" << endl;
         a = rnd.Rannyu(1,ng);
         b = rnd.Rannyu(1,ng);
         c = rnd.Rannyu(1,ng*0.6);
         chr_appo[i].MPermut(a,b,c);
         //cout << a <<" "<< b <<" "<< c << endl;
         //chr_appo[i].Print();
      }

      prob = rnd.Rannyu();
      if(prob < mp[5]){
         //cout << "--CR--" << endl;
         a = rnd.Rannyu(1,ng); // pos
         b = rnd.Rannyu(1,ng); // len
         c = rnd.Rannyu(1,Npop); // altro individuo
         //chr_appo[i].Print();
         chr_appo[i].Crossover(a,b,chr_appo[c]);
         //chr_appo[i].Print();
      }

      // controllo e salvo per la generazione successiva
      chr_appo[i].Check();
      this->Chr[i] = chr_appo[i];
   }
   //delete chr_appo;
}

//========================================================//
//       CHROMOSOME                                       //
//========================================================//

Chromosome :: Chromosome(){} // costruttore
Chromosome :: ~Chromosome(){} // distruttore

Chromosome& Chromosome::operator= (const Chromosome& chr)
{
    // COPIES
    // vettore di geni
    for(int i = 0; i<Ng; i++) Gen[i] = chr.Gen[i];
    // fitness
    Ftn = chr.Ftn;
    //lughezza
    Len = chr.Len;

    // return the existing object so we can chain this operator
    return *this;
}

// methods

// get e set

void Chromosome :: SetGen(int vec[Ng]){
   for(int i = 0; i< Ng; i++) Gen[i] = vec[i];
}

void Chromosome :: SetFitness(double f){ 
   Ftn = f;
}

void Chromosome :: SetLen(double l){ 
   Len = l;
}

double Chromosome :: GetFtn(){
   return Ftn;
}

double Chromosome :: GetLen(){ 
   return Len;
}

int Chromosome :: GetNg(){
   return Ng;
}

int Chromosome :: GetGen(int i){
   return Gen[i];
}


// altri metodi

/// stampa i geni del cromosoma
void Chromosome :: Print(){

   for(int i =0; i<Ng; i++) cout << Gen[i] << " ";
   cout << endl;
}

/// riempie il cromosoma con geni casuali
void Chromosome :: Fill(Random &rnd){
   //if(Gen[0] != 0) {cerr << "Can't fill full chromosome!" << endl; abort();}
   if(Gen[0] != 0) {cerr << "Can't fill full chromosome: I clean it." << endl; this->Empty();}
   
   int r = (int)(rnd.Rannyu()*Ng);
   Gen[0] = 1;
   // prendo interi in ordine
   for(int i = 2; i<Ng+1; i++){
      
      while(Gen[r]!=0){ // finchè la posizione è vuota ...
         r = (int)(rnd.Rannyu()*Ng);
      }
      Gen[r] = i; // ...la riempio...
   }
}

/// svuota il cromosoma
void Chromosome :: Empty(){
   for(int i = 0; i<Ng; i++){
      Gen[i] = 0; 
   }
}


// accessori alle mutazioni

/// Applica pbc quando necessario durante le mutazioni
int Chromosome :: Pbc(int pos){
   if(pos>(Ng-1) && pos%(Ng-1)!=0)    pos = pos%(Ng-1); // non Ng, perchè devo saltare la posizione 0
   if(pos>(Ng-1) && pos%(Ng-1)==0)    pos = Ng-1;
   if(pos < 0)       pos = Ng + pos%Ng;
   return pos;
}

/// Verifica indirettamente (può fallire ma è molto improbabile)
/// che il cromosoma contenga tutti i geni e il primo sia fisso.
/// Verifica che la somma totale dei geni sia corretta.
bool Chromosome :: Check(){

   if(Gen[0] != 1){ // controllo che il primo gene sia sempre 1
      cerr << "Check: something went wrong. First gene is no more 1." << endl;
      exit(1);
   }

   int sum = 0;
   for(int i = 0; i<Ng; i++){
      sum += Gen[i];
   }

   if(Ng == 34 and sum != 595){
      cerr << "Check("+to_string(34)+"): something went wrong. Genes were not conserved through mutations." << endl;
      // << modificare 594 se si usa Ng != 34
      exit(2);
   }

   if(Ng == 50 and sum != 1275){
      cerr << "Check("+to_string(50)+"): something went wrong. Genes were not conserved through mutations." << endl;
      // << modificare 594 se si usa Ng != 34
      exit(2);
   }

   if(Ng == 10 and sum != 55){
      cerr << "Check("+to_string(10)+"): something went wrong. Genes were not conserved through mutations." << endl;
      // << modificare 594 se si usa Ng != 10
      exit(2);
   }
   if(Ng == 20 and sum != 210){
      cerr << "Check("+to_string(20)+"): something went wrong. Genes were not conserved through mutations." << endl;
      // << modificare 594 se si usa Ng != 10
      exit(2);
   }



   return 1;
}


// Mutazioni

/// 1) permuta due geni
void Chromosome :: PairPermut(int pos1, int pos2){
   if(pos1==0) {
      // cerr << "Chromosome::PairPermut : can't move first gene: I move the second one." << endl;
      pos1 = 1;
   }
   if(pos2==0) {
      // cerr << "Chromosome::PairPermut : can't move first gene: I move the second one." << endl;
      pos2 = 1;
   }

   int g1 = Gen[Pbc(pos1)];
   int g2 = Gen[Pbc(pos2)];
   
   Gen[Pbc(pos1)] = g2;
   Gen[Pbc(pos2)] = g1;
}

/// 2) in posizione POS inverte M geni 
void Chromosome :: Inversion(int pos, int m){
   if(pos==0) {
      // cerr << "Chromosome::Inversion : can't move first gene: I move the second one." << endl;
      pos = 1;
   }
   if(pos>(Ng-1)) {
      // cerr << "Chromosome::Inversion : position out of Chromosome dimension. I apply PBC." << endl;
      pos = Pbc(pos);
   }
   if(m > Ng-1) {
      // cerr << "Chromosome::Inversion : too many genes to swap. I set m = 2." << endl; 
      m = 2;
   }

   int block[m];

   // metto da parte il blocco da swappare
   for(int j = 0; j < m; j++)
   {
      block[j] = Gen[Pbc(pos+j)];
   }

   // riempio al contrario il buco lasciato dal blocco
   for(int i = 0; i<m; i++)
   {
      Gen[Pbc(pos+i)] = block[(m-1)-i];
   }

   //delete[] block;
}

/// 3) a partire da una posizione POS, sposta M geni adiacenti in avanti di N posizioni, 
///    (eccetto il primo gene e col vincolo m < Ng-1)
void Chromosome :: Shift(int pos, int m, int n){

   if(pos==0) {
      cerr << "Chromosome::Shift : can't move first gene: I move the second one." << endl;
      pos = 1;
   }
   if(pos>(Ng-1)) {
      cerr << "Chromosome::Shift : position out of Chromosome dimension. I apply PBC." << endl;
      pos = Pbc(pos);
   }
   //if(n<0) {
   //   cerr << "Chromosome::Shift : negative shift. I apply PBC." << endl; 
   //   n = Ng+n%Ng;
   //}
   if(m >= Ng-1) {
      cerr << "Chromosome::Shift : too many genes to move. I set m = 1." << endl; 
      m = 1;
   }

   int block[m];

   // metto da parte il blocco da shiftare
   for(int j = 0; j < m; j++)
   {
      block[j] = Gen[Pbc(pos+j)];
   }

   // shifto n volte i geni complementari per riempire il buco...
   for(int i = 0; i<n; i++)
   {
      Gen[Pbc(pos+i)] = Gen[Pbc(pos+m+i)];
   }

   // ... e rimetto il blocco
   for(int i = 0; i<m; i++)
   {
      Gen[Pbc(pos+n+i)] = block[i];
   }

   //delete[] block;
}

/// 4) shifta tutto di "shift" posizioni
void Chromosome :: Shift2(int shift){
   
   int block[Ng];

   // metto da parte il blocco da shiftare
   for(int j = 1; j < Ng; j++)
   {
      block[j-1] = Gen[j];
   }

   // sovrascrivo il blocco shiftato
   for(int i = 1; i < Ng; i++)
   {
      Gen[Pbc(i+shift)] = block[i-1];
   }
}

/// 5) scambia M geni in posizione POS1 con altrettanti in POS2
void Chromosome :: MPermut(int pos1, int pos2, int m){
   if(pos1==0) {
      cerr << "Chromosome::MPermut : can't move first gene: I move the second one." << endl;
      pos1 = 1;
   }
   if(pos2==0) {
      cerr << "Chromosome::MPermut : can't move first gene: I move the second one." << endl;
      pos2 = 1;
   }
   if(pos1>(Ng-1)) {
      cerr << "Chromosome::MPermut : position out of Chromosome dimension. I apply PBC." << endl;
      pos1 = Pbc(pos1);
   }
   if(pos2>(Ng-1)) {
      cerr << "Chromosome::MPermut : position out of Chromosome dimension. I apply PBC." << endl;
      pos2 = Pbc(pos2);
   }
   if(m > Ng/2) {
      //cerr << "Chromosome::MPermut : too many genes to swap. I set m = 1." << endl; 
      m = 1;
   }
   if(m > abs(pos2-pos1)){ //<< casino
      //cerr << "Chromosome::MPermut : m > abs(pos1-pos2). I set m = abs(pos1-pos2)." << endl; 
      m = abs(pos1-pos2);
   }
   if(pos1>pos2){
      int appo = pos2;
      pos2 = pos1;
      pos1 = appo;
   }

   int block[m];

   // metto da parte il blocco da swappare
   for(int j = 0; j < m; j++)
   {
      block[j] = Gen[Pbc(pos1+j)];
   } 

   // riempio il buco lasciato dal blocco...
   for(int i = 0; i<m; i++)
   {
      Gen[Pbc(pos1+i)] = Gen[Pbc(pos2+i)];
   }

   // ... e rimetto il blocco in pos2
   for(int i = 0; i<m; i++)
   {
      Gen[Pbc(pos2+i)] = block[i];
   }

   //delete[] block;
}

/// 6) fa il crossover di len geni a partire dalla posizione pos, con un cromosoma parent2
void Chromosome :: Crossover(int pos, int len, Chromosome parent2){

   if(pos==0) {
      cerr << "Chromosome::Crossover : can't move first gene: I move the second one." << endl;
      pos = 1;
   }

   int block[len];
   //int block2[len];

   // 1. cut their paths at the same position:
   //    metto da parte i blocchi da swappare
   //    (conservando la prima parte)
   for(int j = 0; j < len; j++)
   {
      block[j] = Gen[Pbc(pos+j)];
      //block2[j] = parent2.Gen[Pbc(pos+j)];
   } 

   // 2. complete the paths with the missing cities adding them in the **order** 
   //    in which they appear in the consort (vale anche se sono separati!):
   int count1 = 0;
   //int count2 = 0;
   
   for(int j = 1; j < Ng; j++)
   {
      // ciclo sui geni nei blocchi 
      for(int i = 0; i < len; i++){
         // quando trovo nell'altro genitore il gene presente
         // nel blocco, lo salvo nel buco scavato all'inizio
         if(parent2.Gen[Pbc(j)] == block[i] ){
            Gen[Pbc(pos+count1)] = parent2.Gen[Pbc(j)];
            count1++;
         }
         //if(Gen[Pbc(j)] == block2[i] ){ // se il crossover modifica solo il cromosoma parent1, questo ciclo non serve
         //   parent2.Gen[Pbc(pos+count2)] = Gen[Pbc(j)];
         //   count2++;
         //}
      }
   } 
}





//========================================================//
//       PROBLEMSET                                       //
//========================================================//

/// Costruttore
Problemset :: Problemset(){}
/// Distruttore
Problemset :: ~Problemset(){}

/// genera coordinate di città in su una circonferenza unitaria
void Problemset :: GenCircCities(){
   
   Random rnd;
   rnd.SetSeed();

   double phi;

   for(int i = 0; i<Ncit; i++){
      phi = rnd.Rannyu()*2*M_PI;
      Xcit[i] = cos(phi);
      Ycit[i] = sin(phi);
   }
}

/// genera coordinate di città in un quadrato di lato unitario
void Problemset :: GenSquareCities(){
   
   Random rnd;
   rnd.SetSeed();

   double pos;

   for(int i = 0; i<Ncit; i++){
      pos = rnd.Rannyu()*2-1;
      Xcit[i] = pos;
      pos = rnd.Rannyu()*2-1;
      Ycit[i] = pos;
   }
}

/// stampa le coordinate delle città nell'ordine indicato da un cromosoma
void Problemset :: PrintCities(int generation, Chromosome chr){

   ofstream stream;
   stream.open("cit/"+to_string(generation)+"citycoord.out");

   int seq;
   for(int i = 0; i<Ncit; i++){
      // stampo 
      seq = chr.GetGen(i)-1; // nell'ordine indicato dal cromosoma
      stream << Xcit[seq]  << " " << Ycit[seq] << endl;
   }
   stream.close();
}

/// stampa la lunghezza media della migliore metà di individui in una popolazione
void Problemset :: PrintBestsLenAve(int generation, int part, Population pop){
   
   ofstream str;
   str.open("bestLen/BLAv.out", std::ios_base::app);

   // ho una variabile BestsLenAve in pop che tiene 
   // conto della media nella metà migliore 
   // (viene riempita quando faccio le mutazioni)

   str << generation << " " << pop.BestsLenAve << endl;
   str.close();
}

/// Stampa la lunghezza del miglior individuo di una popolazione
void Problemset :: PrintBestLen(int generation, Population pop){

   ofstream str;
   str.open("bestLen/BL.out", std::ios_base::app);

   // ho una variabile BestLen in pop che tiene 
   // conto della lunghezza migliore 
   // (viene riempita quando faccio le mutazioni)

   str << generation << " " << pop.BestLen << endl;
   str.close();
}

/// prende un cromosoma e con l'ordine dei geni calcola la distanza tra le città e la fitness
void Problemset :: EvalFitness(Chromosome &chr){
   double fitn = 0;
   double accu = 0;
   double acculen = 0;
   double len2 = 0;
   int ng = chr.GetNg();

   // segmenti dalla prima all'ultima
   for(int i = 0; i< ng-1; i++){

      len2 = pow( Xcit[chr.GetGen(i+1)-1] -Xcit[chr.GetGen(i)-1] ,2) 
           + pow( Ycit[chr.GetGen(i+1)-1] -Ycit[chr.GetGen(i)-1] ,2);
      accu += len2;
      acculen += sqrt(len2);
   }
   // distanza tra la prima e l'ultima
   len2 = pow(Xcit[chr.GetGen(0)-1]-Xcit[chr.GetGen(ng-1)-1],2) 
         + pow(Ycit[chr.GetGen(0)-1]-Ycit[chr.GetGen(ng-1)-1],2);
   accu += len2;
   acculen += sqrt(len2);

   fitn = 1.0/accu;
   chr.Ftn = fitn;
   chr.Len= acculen;

}

/// prende una popolazione e applica EvalFitness su tutti i suoi individui
void Problemset :: EvalAll(Population &pop){ 
   for(int i = 0; i < pop.GetNPop(); i++){
      EvalFitness(pop.Chr[i]);
   }
}

/// prende una popolazione e mette in ordine i suoi individui dal peggiore al migliore
void Problemset :: SortPop(Population *pop){
   quickSort(pop, 0, pop->GetNPop()-1);
}


//========================================================//
//       OTHER FUNCTIONS                                  //
//========================================================//

/// funzione necessaria in quicksort
int partition(Population *pop, int low, int high){
    
   // prendo la ftn dell'elemento più a destra
   double pivot = pop->Chr[high].Ftn;

   // puntatore all'elemento più grande (ipotizzo sia a sx) 
   int i = (low-1);

   // comparo ogni ftn del vettore con la ftn pivot
   for (int j = low; j < high; j++) {
      if (pop->Chr[j].Ftn <= pivot) {
         // se ftn è minore, scambio l'intero chr con quello puntato da i
         i++;

         // swap
         Chromosome t = pop->Chr[i];
         pop->Chr[i] = pop->Chr[j];
         pop->Chr[j] = t;
      }
   }

   // scambio pivot con l'elemento puntato da i
   Chromosome tt = pop->Chr[i+1];
   pop->Chr[i+1] = pop->Chr[high];
   pop->Chr[high] = tt;

   // return the partition point
   return (i + 1);
}

/// funzione necessaria in SortPop
void quickSort(Population *pop, int low, int high) {
  if (low < high) {

    int pi = partition(pop, low, high);

    quickSort(pop, low, pi - 1); 
    quickSort(pop, pi + 1, high); 
   }
}
