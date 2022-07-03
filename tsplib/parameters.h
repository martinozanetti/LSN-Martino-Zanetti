
static const int NindPop = 3000;      // number of idividuals in population [good:3000]
static const int Ngenes = 50;        // number of genes in a chromosome = number of cities in the problemset
static const int NGeneration = 500; // numero di generazioni consecutive
static const double pmut[6] = {
    0.1, // PairPermut
    0.1, // Invertion
    0.2, // Shift
    0.2, // Shift2
    0.4, // MPermut
    0.6  // Crossover
    };
static const double expon = 6.0;