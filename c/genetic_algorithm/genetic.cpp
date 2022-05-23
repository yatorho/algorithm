// Copyright by 2022.3 yatorho
// author yatorho

#include "genetic.h"
#include <cstdlib>
#include <memory.h>
#include <random>

namespace genetic {

template<typename Dtype>
Dtype *GeneticAlgorithm<Dtype>::start() {
  Dtype *fitness;
  Dtype *genes_list;
  Dtype *genes_next;

  fitness = static_cast<Dtype *>(malloc(_po_size * sizeof(Dtype)));
  genes_list = static_cast<Dtype *>(
          malloc(_po_size * _geneLength * sizeof(Dtype)));
  genes_next = static_cast<Dtype *>(
          malloc(_po_size * _geneLength * sizeof(Dtype)));

  _init_population(fitness, genes_list);


  for (uint32_t ev = 0; ev < _ev_algebra; ev++) {

    _sort_by_fitness(fitness, genes_list);

    std::cout << "Best:" << fitness[_po_size - 1] << std::endl;
    std::cout << "Average: " << _average(fitness) << std::endl;
    std::cout << "Times: " << ev << std::endl;

    _roulette_select(fitness, genes_list, genes_next, 0.2, 0.05);

    _crossing_operator(genes_next);

      _mutating_operator(genes_next);

      memcpy(genes_list, genes_next, _po_size * _geneLength * sizeof(Dtype));

      _fitting_operator(genes_list, fitness);
  }

  Dtype *best_genes;
  best_genes = static_cast<Dtype *>(malloc(_geneLength * sizeof(Dtype)));

  memcpy(best_genes, genes_list + (_po_size - 1) * _geneLength,
         _geneLength * sizeof(Dtype));

  free(fitness);
  free(genes_list);
  free(genes_next);


  return best_genes;
}


template<typename Dtype>
GeneticAlgorithm<Dtype>::GeneticAlgorithm(uint32_t Po_Size, uint32_t EV_Algebra,
                                          float P_crossover, float P_mutation,
                                          uint32_t gene_len)
        : _po_size(Po_Size),
          _ev_algebra(EV_Algebra),
          _p_crossover(P_crossover),
          _p_mutation(P_mutation),
          _geneLength(gene_len) {
  gen = std::mt19937(rd());
  dis = std::uniform_real_distribution<double>(0, 1);
  MessageBox mb;
}

template<typename Dtype>
GeneticAlgorithm<Dtype>::GeneticAlgorithm(uint32_t gene_len)
        : _po_size(100ul),
          _ev_algebra(100ul),
          _p_crossover(0.7f),
          _p_mutation(0.05f),
          _geneLength(gene_len) {
  gen = std::mt19937(rd());
  dis = std::uniform_real_distribution<double>(0, 1);
  MessageBox mb;
}

template<typename Dtype>
void GeneticAlgorithm<Dtype>::_sort_by_fitness(Dtype *fitness, Dtype *genes) {
  Dtype temp_f;
  Dtype *temps_gs;
  temps_gs = static_cast<Dtype *>(malloc(_geneLength * sizeof(Dtype)));

  for (uint32_t i = 0; i < _po_size - 1; i++) {
    for (uint32_t j = i + 1; j < _po_size; j++) {
      if (fitness[i] > fitness[j]) {
        temp_f = fitness[j];
        fitness[j] = fitness[i];
        fitness[i] = temp_f;

        memcpy(temps_gs, genes + offset(j, 0), _geneLength * sizeof(Dtype));
        memcpy(genes + offset(j, 0), genes + offset(i, 0),
               _geneLength * sizeof(Dtype));
        memcpy(genes + offset(i, 0), temps_gs, _geneLength * sizeof(Dtype));
      }
    }
  }
  free(temps_gs);
}

template<typename Dtype>
uint32_t
GeneticAlgorithm<Dtype>::offset(uint32_t index_0, uint32_t index_1) const {
  return index_0 * _geneLength + index_1;
}

template<typename Dtype>
Dtype *GeneticAlgorithm<Dtype>::start(Dtype (*func)(Dtype *, uint32_t),
                                      Dtype *(*get_genes)()) {
  return nullptr;
}

template<typename Dtype>
Dtype *GeneticAlgorithm<Dtype>::_start() {
  return nullptr;
}

template<typename Dtype>
void
GeneticAlgorithm<Dtype>::_init_population(Dtype *fitness, Dtype *genes_list) {
  for (uint32_t p = 0; p < _po_size; p++) {
    EncodeGenes(genes_list + p * _geneLength);
    fitness[p] = function(genes_list + p * _geneLength);
  }
}

template<typename Dtype>
void GeneticAlgorithm<Dtype>::sort_by_fitness(Dtype *fitness, Dtype *genes) {
  _sort_by_fitness(fitness, genes);
}

template<typename Dtype>
Dtype GeneticAlgorithm<Dtype>::_average(const Dtype *fitness) {
  auto sum = (Dtype) 0.0;
  for (uint32_t i = 0; i < _po_size; i++) {
    sum += fitness[i];
  }
  return sum / _po_size;
}

template<typename Dtype>
void GeneticAlgorithm<Dtype>::EncodeGenes(Dtype *genes) {}

template<typename Dtype>
void
GeneticAlgorithm<Dtype>::_roulette_select(const Dtype *fitness,
                                          const Dtype *genes_list,
                                          Dtype *genes_next,
                                          double weed_rate,
                                          double saved_head) {
  saved_head = 1.0 - saved_head;

  double roulette = 0.0;
  for (auto i = (uint32_t) (weed_rate * _po_size); i < _po_size; i++) {
    roulette += exp((double) fitness[i]);
  }

  uint32_t len_flag = 0ul;
  for (auto i = (uint32_t) (weed_rate * _po_size); i < _po_size; i++) {
    double prob = exp((double) fitness[i]) / roulette;
    auto num = (uint32_t) (round((uint32_t) (_po_size * saved_head) * prob));
    for (uint32_t j = 0; j < num; j++) {
      if (len_flag < (uint32_t) (_po_size * saved_head)) {
        memcpy(genes_next + len_flag * _geneLength,
               genes_list + i * _geneLength, _geneLength * sizeof(Dtype));
        len_flag += 1ul;
      }
    }
  }

  for (uint32_t i = len_flag; i < _po_size; i++) {
    memcpy(genes_next + len_flag * _geneLength, genes_list + i * _geneLength,
           _geneLength * sizeof(Dtype));
    len_flag += 1ul;
  }

}

template<typename Dtype>
void GeneticAlgorithm<Dtype>::_crossing_operator(Dtype *genes_list) {

  uint32_t crosser_mark = _geneLength / 2;
  for (uint32_t i = 0; i < _po_size - 1; i++) {
    if (dis(gen) < _p_crossover) {
      int which_gene = static_cast<int>(round(dis(gen) * (_po_size - 2)));

      if (dis(gen) < 0.5) {
        memcpy(genes_list + which_gene * _geneLength,
               genes_list + i * _geneLength, crosser_mark * sizeof(Dtype));
      } else {
        memcpy(genes_list + which_gene * _geneLength + crosser_mark,
               genes_list + i * _geneLength + crosser_mark,
               crosser_mark * sizeof(Dtype));
      }
    }
  }
}

template<typename Dtype>
void GeneticAlgorithm<Dtype>::_mutating_operator(Dtype *genes_list) {

  auto *mutation_gens = static_cast<Dtype *>(malloc(
          _geneLength * sizeof(Dtype)));

  for (uint32_t p = 0; p < _po_size - 1; p++) {
    EncodeGenes(mutation_gens);
    for (uint32_t g = 0; g < _geneLength; g++) {
      if (dis(gen) < _p_mutation) {
        genes_list[p * _geneLength + g] = mutation_gens[g];
      }
    }
  }
  free(mutation_gens);
}

template<typename Dtype>
void
GeneticAlgorithm<Dtype>::_fitting_operator(const Dtype *genes, Dtype *fitness) {
  for (uint32_t i = 0; i < _po_size; i++) {
    fitness[i] = function(genes + i * _geneLength);
  }
}

template<typename Dtype>
Dtype GeneticAlgorithm<Dtype>::function(const Dtype *genes) {
  return static_cast<Dtype>(0);
}

template<typename Dtype>
Dtype GeneticAlgorithm<Dtype>::max(const Dtype *arr, uint32_t len) {
  Dtype max_value = 0.0;
  for (uint32_t i = 0; i < len; i++) {
    if (max_value < arr[i]) {
      max_value = arr[i];
    }
  }
  return max_value;
}

template<typename Dtype>
int GeneticAlgorithm<Dtype>::index_value(const Dtype *arr, Dtype value,
                                         uint32_t len) {
  for (uint32_t i = 0; i < len; i++) {
    if (arr[i] == value) {
      return i;
    }
  }
  return -1;
}


template<typename Dtype>
GeneticAlgorithm<Dtype>::GeneticAlgorithm() = default;

template float GeneticAlgorithm<float>::function(const float *genes);

template double GeneticAlgorithm<double>::function(const double *genes);

template void
GeneticAlgorithm<float>::sort_by_fitness(float *fitness, float *genes);

template void
GeneticAlgorithm<double>::sort_by_fitness(double *fitness, double *genes);

template float *GeneticAlgorithm<float>::start();

template double *GeneticAlgorithm<double>::start();


template GeneticAlgorithm<double>::GeneticAlgorithm(uint32_t gene_len);

template GeneticAlgorithm<float>::GeneticAlgorithm(uint32_t gene_len);

template GeneticAlgorithm<int>::GeneticAlgorithm(uint32_t gene_len);

template GeneticAlgorithm<uint32_t>::GeneticAlgorithm(uint32_t gene_len);

template GeneticAlgorithm<double>::GeneticAlgorithm(uint32_t Po_Size,
                                                    uint32_t EV_Algebra,
                                                    float P_crossover,
                                                    float P_mutation,
                                                    uint32_t gene_len);

template GeneticAlgorithm<float>::GeneticAlgorithm(uint32_t Po_Size,
                                                   uint32_t EV_Algebra,
                                                   float P_crossover,
                                                   float P_mutation,
                                                   uint32_t gene_len);

template GeneticAlgorithm<int>::GeneticAlgorithm(uint32_t Po_Size,
                                                 uint32_t EV_Algebra,
                                                 float P_crossover,
                                                 float P_mutation,
                                                 uint32_t gene_len);

template GeneticAlgorithm<uint32_t>::GeneticAlgorithm(uint32_t Po_Size,
                                                      uint32_t EV_Algebra,
                                                      float P_crossover,
                                                      float P_mutation,
                                                      uint32_t gene_len);

template void
GeneticAlgorithm<double>::_sort_by_fitness(double *fintness, double *genes);

template void
GeneticAlgorithm<float>::_sort_by_fitness(float *fintness, float *genes);

template void
GeneticAlgorithm<float>::_init_population(float *fitness, float *genes_list);

template void
GeneticAlgorithm<double>::_init_population(double *fitness, double *genes_list);

template float GeneticAlgorithm<float>::_average(const float *fitness);

template double GeneticAlgorithm<double>::_average(const double *fitness);

template void GeneticAlgorithm<float>::_roulette_select(const float *fitness,
                                                        const float *genes_list,
                                                        float *genes_next,
                                                        double weed_rate,
                                                        double saved_head);

template void GeneticAlgorithm<double>::_roulette_select(const double *fitness,
                                                         const double *genes_list,
                                                         double *genes_next,
                                                         double weed_rate,
                                                         double saved_head);

template void GeneticAlgorithm<float>::_crossing_operator(float *genes_next);

template void GeneticAlgorithm<double>::_crossing_operator(double *genes_next);

template void GeneticAlgorithm<float>::_mutating_operator(float *genes_list);

template void GeneticAlgorithm<double>::_mutating_operator(double *genes_list);

template void
GeneticAlgorithm<float>::_fitting_operator(const float *genes, float *fitness);

template void GeneticAlgorithm<double>::_fitting_operator(const double *genes,
                                                          double *fitness);

template float GeneticAlgorithm<float>::max(const float *arr, uint32_t len);

template double GeneticAlgorithm<double>::max(const double *arr, uint32_t len);

template int GeneticAlgorithm<float>::index_value(const float *arr, float value,
                                                  uint32_t len);

template int
GeneticAlgorithm<double>::index_value(const double *arr, double value,
                                      uint32_t len);

}  // namespace genetic
