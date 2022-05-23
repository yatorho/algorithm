#include "genetic.h"
#include <random>
#include <ctime>

// #include <cstdlib>
// #include <ctime>


static uint32_t GENES_LENGTH = 10;
static uint32_t PO_SIZE = 1000;
static uint32_t EV_ALGEBRA = 200;
static float P_CROSSOVER = 0.7f;
static float P_MUTATION = 0.1f;

class MyGenetic : public genetic::GeneticAlgorithm<double> {
public:
  explicit MyGenetic(uint32_t gen_length) : genetic::GeneticAlgorithm<double>(
          gen_length) {
    srand((int) time(nullptr));
  }

  MyGenetic(uint32_t po_size, uint32_t ev_algebra, float p_crossover,
            float p_mutation,
            uint32_t gen_length) : genetic::GeneticAlgorithm<double>(po_size,
                                                                     ev_algebra,
                                                                     p_crossover,
                                                                     p_mutation,
                                                                     gen_length) {
    srand((int) time(nullptr));
  }

  double function(const double *genes) override {
    double res = 0.0;
    for (uint32_t i = 0; i < GENES_LENGTH; i++) {
      res += (genes[i] - i) * (genes[i] - i);
    }
    return -res;
  }


  void EncodeGenes(double *genes) override {

    for (uint32_t i = 0; i < GENES_LENGTH; i++) {
      double n = (rand() / (double) (RAND_MAX) - 0.5) * 30;
      // std::cout << n << std::endl;
      genes[i] = n;
    }
  }
};


int main() {
  double *best;

  MyGenetic ga(PO_SIZE, EV_ALGEBRA, P_CROSSOVER, P_MUTATION, GENES_LENGTH);
  best = ga.start();

  for (uint32_t i = 0; i < GENES_LENGTH; i++) {
    std::cout << best[i] << " ";
  }
  std::cout << std::endl;

  free(best);
  return 0;
}