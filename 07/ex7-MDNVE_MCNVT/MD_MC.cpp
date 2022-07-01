/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ex7

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

string phase = "no_phase";

int main(int argc, char* argv[])
{ 
  
  API((int)argc, argv);
  Input();

  //int nconf = 1;
  // ciclo sui blocchi
  for(int iblk=1; iblk <= nblk; iblk++) 
  {
    Reset(iblk); 
    //ciclo nel blocco
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure(); // misurato con correzioni di coda
      Accumulate(); // Update block averages
      /*
      if(istep%10 == 0){
      //  ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! cioè ti riempie facilmente la RAM!
        nconf += 1;
      }*/
      if (verbose and istep%1000==0){
        cout << "Block: " << iblk << "/"<< nblk<< ", step: " << istep << endl;
      }
    }
    Averages(iblk);   //Print results for current block
  }

  // NON registro la configuazione finale per non rischiare di dover riequilibrare
  ConfFinal(); //Write final configuration  

  return 0;
}

// ===================== FUNZIONI ====================================

void API(int argc, char *argv[]){
  
  if (argc == 1)
  { int tentativi = 0;
    while (phase != "solid" && phase != "liquid" && phase != "gas")
    {
      if (tentativi > 0) cout << endl << "Error: \'" << phase << "\' is not a valid phase" << endl;
      cout << "Please type a phase (or directly run \"./NVE_NVT.exe <phase>\"): solid, liquid or gas? ";
      cin >> phase;
      tentativi ++;
    }
  }
  else 
  { 
    phase = (string)argv[1];
    while (phase != "solid" && phase != "liquid" && phase != "gas")
    { 
      cout << endl << "Error: \'" << phase << "\' is not a valid phase" << endl;
      cout << "Usage: ./NVE_NVT.exe <phase>, whith phase={solid, liquid, gas}." << endl;
      cout << "Please retype a phase (or directly run \"./NVE_NVT.exe <phase>\"): solid, liquid or gas? ";
      cin >> phase;
    }
  }

  cout << endl << "Did you remember cleaning .dat files (./clean"+phase+".sh)? [y/n]" << endl;
  char resp;
  while (resp != 'y' && resp!= 'n')
  {
    cin >> resp;
    if (resp == 'y' ) { }
    else if (resp == 'n') { 
      cout << "Allora ciao!" << endl;
      exit(0); 
    }
  }
}

// vvvvv nessuna modifica per la g(r)
void Input(void) 
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << endl;
  cout << ",-------------------------------,    " << endl;
  cout << "| Classic Lennard-Jones fluid   |    " << endl;
  cout << "| MD(NVE) / MC(NVT) simulation  |    " << endl;
  cout << "'-------------------------------'    " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl;
  cout << "The program uses Lennard-Jones units " << endl << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("../../random/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input."+ phase);

  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open("seed.out");
  else Seed.open("../../random/seed.in");
  
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;

  ReadInput >> npart;
  ReadInput >> rho;
  vol = (double)npart/rho; // calcolo volume a partire dalla densità
  box = pow(vol,1.0/3.0);
  ReadInput >> rcut;
  ReadInput >> delta; // timestep in unità naturali: dt*
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> verbose;

  cout << ",------------- Input parameters --------------," << endl;
  cout << "| Temperature                         = " << temp << endl;
  cout << "| Number of particles                 = " << npart << endl;
  cout << "| Density of particles                = " << rho << endl;
  cout << "| Volume of the simulation box        = " << vol << endl;
  cout << "| Edge of the simulation box          = " << box << endl;
  cout << "| Cutoff of the interatomic potential = " << rcut << endl;
  cout << "'---------------------------------------------'" << endl << endl;
    
  cout << "The program performs Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter =              " << delta << endl;
  cout << "Number of blocks =             " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements: per ora ho 4 indici
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  iw = 4; //Pressure: non uso ip perchè si scontra con altre variabili locali (es: in Boltzmann, Force...)
  n_props = 5; //Number of observables

//Read initial configuration
  if(restart) // se devo ripartire da dove ero rimasto, carico l'ultima configurazione
  {
    if(iNVET){
      ReadConf.open("init/" + phase + "/config.out");
      ReadVelocity.open("init/" + phase + "/velocity.out");
      for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
    }else if (!iNVET){
      ReadConf.open("init/MD/" + phase + "/config.out");
      ReadVelocity.open("init/MD/" + phase + "/velocity.out");
      for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
    }
  }
  else // altrimenti apro la config cristallo perfetto e preparo le velocità, compatibilmente con la simul che devo fare
  {
    ReadConf.open("init/config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0}; // velocity sum
    for (int i=0; i<npart; ++i) // creo le velocità, con distr. gaussiana (maxwell-boltzmann)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
      // problema: così gneero una (piccola) velocità di drift, che ci portiamo avanti per sempre come en cinetica in più non realistica ...
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart; //... accumulo le velocità di drift
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0]; //... e le sottraggo
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]; // ... allora calcolo un fattore di riscalamento per "contro-correggere" il conto della temperatura
    }
    sumv2 /= (double)npart; // normalizzo
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs; // riscalo: a questo punto è garantito che non ho il drift e le velocità sono corrette, per ogni temperatura data
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i) //settaggio delle posizioni di tutte le particelle
  {
    ReadConf >> x[i] >> y[i] >> z[i]; // è input.in-out a seconda che stia ripartendo da zero o da un passo precedente
    x[i] = Pbc( x[i] * box ); // NOTA: LE COORDINATE SONO SCRITTE IN UNITÀ DEL LATO: DA -0,5 A +0,5, USATE NEL FILE CONFIG: 
    y[i] = Pbc( y[i] * box ); // ... quindi la priuma cosa che faccio è moltiplicarlo per il lato
    z[i] = Pbc( z[i] * box ); // per scrupolo, infilo subito il valore trovato nella funzione delle PBC
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i) // applico algoritmo (verlet o MC)
  {
    if(iNVET) // se voglio monte carlo. Nota: int invet = 0 corrisponde a false, int invet = 1 corrisponde a true
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else // se voglio MD simulation. 
    {
      xold[i] = Pbc(x[i] - vx[i] * delta); // alg. di verlet: compute old config from velocity
      yold[i] = Pbc(y[i] - vy[i] * delta); // ATTENZIONE: questa traiettoria non è puro verlet, perchè mischia velocità e posizione. 
      zold[i] = Pbc(z[i] - vz[i] * delta); // ... quindi è un po meno preciso, ma per quel che facciamo noi va benissimo.
                                           // ... se uno vuole può aggiustare. ma non è così rilevante
    }
  }
  
//Calculate tail corrections <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  v_tail = 8 * pi * rho * (1.0/(9*pow(rcut,9))-1.0/(3*pow(rcut,3))) * npart;
  w_tail = 32* pi * rho * (1.0/(9*pow(rcut,9))-1.0/(6*pow(rcut,3)))*3*npart;

//Evaluate properties of the initial configuration, tanto per vedere le proprietà iniziali
  Measure();

//Print initial values for measured properties
  cout << ",-------- Initial configuration -------," << endl;
  cout << "| Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "| Initial temperature      = " << walker[it] << endl;
  cout << "| Initial pressure         = " << walker[iw] << endl;
  cout << "'--------------------------------------'" << endl;

  return;
}

void Measure() //Properties measurement
{
  // azzero tutto
  double v = 0.0, kin=0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;
  for(int i=0; i<ng; i++) gdr[i] = 0; // azzero l'istogramma 

  // cycle over pairs of particles
  for (int i=0; i<npart-1; ++i) // calculate tot potential energy
  {
    for (int j=i+1; j<npart; ++j)
    {
      // distance i-j in pbc = Min Img Convention
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut) // cutoff
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6); // en potenziale
        wij = 48.0/pow(dr,12) - 24.0/pow(dr,6); // viriale
        v += vij;
        w += wij;
      }

      // riempio l'istogramma della gdr
      min_dist = box/2.0; // risoluzione dell'istogramma della gdr
      bin_index = ng*(dr/min_dist); // posizione nell'istogramma

      if(dr<min_dist){
        gdr[ bin_index ] += 2; //<<<<<<<<<<<< gdr
      }
    }
  }  

  for (int i=0; i<npart; ++i){
    kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); // calculate tot kin energy
  }

// walker è il vettore che tiene tutte le quantità: le registro

  walker[iv] = 4.0 * v + v_tail;           // Potential energy, vedi slide 43 (lez.5?) <<<<< aggiunto tail corrections // da dividere per Npart più avanti
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature: dall'en cinetica
  walker[iw] = rho * walker[it] + (w + w_tail)/(3.0*vol);  // Pressure: dal viriale. <<<<< idem

  return;
}

void Averages(int iblk) //Print results for current block. Alla fine del blocco prende ogni proprietà, fa media ed errori standard.
{
    
    ofstream Epot, /*Ekin, Etot, Temp,*/ Pres, Gdr; 
    
    if (verbose && iblk%1000==0)
    {
      cout << "Block number " << iblk << "/" << nblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
      cout << "----------------------------" << endl << endl;
    }

    Epot.open(phase + "/output_epot.dat",ios::app);
    Pres.open(phase + "/output_pres.dat",ios::app);
    Gdr.open(phase + "/output_gdr.dat", ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = blk_av[iw]/blk_norm; //Pressure: dividere per npart? no, essendo una grandezza intensiva
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

    double gdr_norm=0;              // GDR
    for(int i=0; i<ng; i++)         // normalizzo le medie della gdr
    {
      double r = i*(min_dist)/ng; // distanza tra due particelle
      double dr =  (min_dist)/ng;  // risoluzione della gdr
      gdr_norm=   rho*npart 
                  *4.*M_PI/3.0 
                  *(  pow(r+dr, 3)
                     -pow(  r , 3)
                   );
      gdr_ave[i] = double(gdr_ave[i]/(double)blk_norm)/gdr_norm; // medie su gdr.
    }

//Potential energy per particle
  Epot << iblk <<  " " << stima_pot << " " << glob_av[iv]/(double)iblk << " " << err_pot << endl;
//Pressure
  Pres << iblk <<  " " << stima_pres << " " << glob_av[iw]/(double)iblk << " " << err_press << endl;
//Gdr
  for(int i=0; i<ng; i++){
    Gdr << gdr_ave[i] << " "; // stampo le medie di blocco su ogni bin
  }
  Gdr << endl;

  Epot.close();
  Pres.close();
  Gdr.close();
}

// qui sotto non ci sono dipendenze dall'aggiunta della stima di pressione

void Move() // muove da una configurazione a quella successiva
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );
      
      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted++;

      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted++;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme: ciclo che muove tutte le particelle

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) ); // algoritmo di verlet
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) ); // NOTA: PBC, perchè le particelle potrebbero essere scappate dalla scatola
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta); // misuro le velocità in pbc (rischia di divergere la velocità attraverso una parete)
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i]; // aggiorno la conf vechia con quella nuova
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew; // e salvo quella nuova
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0; // variabile che servono solo per monte carlo, per evitare dei NAN che vengono da dopo
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip) // restituisce l'energia pot di una particella tra le altre con pot di lenn-john (?)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){ // cicla su tutte le particelle tranne quella considerata, che non produce forza su se stessa
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr); // modulo della distanza

      if(dr < rcut){ // applico il cutof
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f; // già in unità naturali
}

/// reset block averages
void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
      for(int i=0; i<n_props; ++i)
      {
        glob_av[i] = 0;
        glob_av2[i] = 0;
      }
      for(int i=0; i<ng; i++) //<<<<<<< gdr
      {
        gdr[i]=0;
        gdr_ave[i]=0;
      }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   for(int i=0; i<ng; i++) //<<<<<<< gdr
   {
     gdr[i]=0;
     gdr_ave[i]=0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
   
}

void Accumulate(void) //Update block averages
{
   for (int i=0; i<100; i++)
   {
    gdr_ave[i]=gdr_ave[i]+gdr[i];
    gdr[i]=0;
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void ConfFinal(void) // salva su file le configurazioni spaziali e di velocità
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open(phase + "/config.out");
  WriteVelocity.open(phase + "/velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
