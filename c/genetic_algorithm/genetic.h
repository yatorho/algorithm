// Copyright by 2022.2.21 trchime
// author: yatorho

#ifndef GENETIC_ALGORITHM_GENETIC_H
#define GENETIC_ALGORITHM_GENETIC_H

#include <cstdint>
#include "message.h"
#include <random>


namespace genetic {

template<typename Dtype>
class GeneticAlgorithm {
public:
  virtual Dtype function(const Dtype *genes);

  virtual void EncodeGenes(Dtype *genes);

  GeneticAlgorithm(uint32_t Po_Size, uint32_t EV_Algebra, float P_crossover,
                   float P_mutation,
                   uint32_t gene_len);

  explicit GeneticAlgorithm(uint32_t gene_len);

  void sort_by_fitness(Dtype *fitness, Dtype *genes);

  Dtype *start();

  Dtype *start(Dtype (*func)(Dtype *, uint32_t), Dtype *(*get_genes)());

  Dtype *_start();

  Dtype max(const Dtype *arr, uint32_t len);

  int index_value(const Dtype *arr, Dtype value, uint32_t len);

private:

  GeneticAlgorithm();

  void _init_population(Dtype *fitness, Dtype *genes_list);

  void _sort_by_fitness(Dtype *fitness, Dtype *genes);

  Dtype _average(const Dtype *fitness);

  void _roulette_select(const Dtype *fitness, const Dtype *genes_list,
                        Dtype *genes_next,
                        double weed_rate, double saved_head);

  void _crossing_operator(Dtype *genes_next);

  void _mutating_operator(Dtype *genes_list);

  void _fitting_operator(const Dtype *genes, Dtype *fitness);


  uint32_t _po_size;
  uint32_t _ev_algebra;
  float _p_crossover;
  float _p_mutation;

  uint32_t _geneLength;

  inline uint32_t offset(uint32_t index_0, uint32_t index_1) const;

protected:
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<double> dis;
};
}
#endif  // GENETIC_ALGORITHM_GENETIC_H
