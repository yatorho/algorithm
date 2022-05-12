#include "particle_v2.h"
#include <iostream>
#include <stdint.h>
#include <thread>
#include <vector>
#include <windows.h>

namespace particle {

uint32_t find_best_position_down(const float *arr, uint32_t po_size) {
  uint32_t index = 0ul;
  float best_value = arr[0];
  float temp_value;
  for (uint32_t i = 1; i < po_size; i++) {
    temp_value = arr[i];
    if (best_value > temp_value) {
      best_value = temp_value;
      index = i;
    }
  }
  return index;
}

uint32_t find_best_position_up(const float *arr, uint32_t po_size) {
  uint32_t index = 0ul;
  float best_value = arr[0];
  float temp_value;
  for (uint32_t i = 1; i < po_size; i++) {
    temp_value = arr[i];
    if (best_value < temp_value) {
      best_value = temp_value;
      index = i;
    }
  }
  return index;
}

void init_particle_down(Particle &pco, float (*op_fun)(const float *),
                        void (*define_solution_domain)(float *),
                        void (*define_velocity_range)(float *)) {
  pco.velocity_vector =
      static_cast<float *>(malloc(pco.po_size * pco.pos_len * sizeof(float)));
  pco.position_vector =
      static_cast<float *>(malloc(pco.po_size * pco.pos_len * sizeof(float)));
  pco.p_best_p =
      static_cast<float *>(malloc(pco.po_size * pco.pos_len * sizeof(float)));
  pco.g_best_p = static_cast<float *>(malloc(pco.pos_len * sizeof(float)));

  pco.p_func_value = static_cast<float *>(malloc(pco.po_size * sizeof(float)));

  for (uint32_t i = 0; i < pco.po_size; i++) {
    define_solution_domain(pco.position_vector + i * pco.pos_len);
    define_velocity_range(pco.velocity_vector + i * pco.pos_len);

    pco.p_func_value[i] = op_fun(pco.position_vector + i * pco.pos_len);
  }

  memcpy(pco.p_best_p, pco.position_vector,
         pco.po_size * pco.pos_len * sizeof(float));
  uint32_t g_best_index =
      find_best_position_down(pco.p_func_value, pco.po_size);
  memcpy(pco.g_best_p, pco.position_vector + g_best_index * pco.pos_len,
         pco.pos_len * sizeof(float));
}

void init_particle_up(Particle &pco, float (*op_fun)(const float *),
                      void (*define_solution_domain)(float *),
                      void (*define_velocity_range)(float *)) {
  pco.velocity_vector =
      static_cast<float *>(malloc(pco.po_size * pco.pos_len * sizeof(float)));
  pco.position_vector =
      static_cast<float *>(malloc(pco.po_size * pco.pos_len * sizeof(float)));
  pco.p_best_p =
      static_cast<float *>(malloc(pco.po_size * pco.pos_len * sizeof(float)));
  pco.g_best_p = static_cast<float *>(malloc(pco.pos_len * sizeof(float)));

  pco.p_func_value = static_cast<float *>(malloc(pco.po_size * sizeof(float)));

  for (uint32_t i = 0; i < pco.po_size; i++) {
    define_solution_domain(pco.position_vector + i * pco.pos_len);
    define_velocity_range(pco.velocity_vector + i * pco.pos_len);

    pco.p_func_value[i] = op_fun(pco.position_vector + i * pco.pos_len);
  }

  memcpy(pco.p_best_p, pco.position_vector,
         pco.po_size * pco.pos_len * sizeof(float));
  uint32_t g_best_index = find_best_position_up(pco.p_func_value, pco.po_size);
  memcpy(pco.g_best_p, pco.position_vector + g_best_index * pco.pos_len,
         pco.pos_len * sizeof(float));
}

void show_array(const float *arr, uint32_t len1, uint32_t len2) {
  for (uint32_t i = 0; i < len1; i++) {
    for (uint32_t j = 0; j < len2; j++) {
      std::cout << arr[i * len2 + j] << " ";
    }
    std::cout << std::endl;
  }
}

void start_down(Particle &pco, float (*op_fun)(const float *)) {
  for (uint32_t e = 0; e < pco.ev_count; e++) {
    update_down(pco, op_fun);
    std::cout << "EV: " << e << std::endl;
    if (pco.LOG_BEST_DATA)
      std::cout << "BEST: " << op_fun(pco.g_best_p) << std::endl;
    if (pco.LOG_AVERAGE_DATA)
      std::cout << "AVERAGE: "
                << average(pco.position_vector, op_fun, pco.po_size,
                           pco.pos_len)
                << std::endl;
    std::cout << "========================" << std::endl;
  }
}
void start_up(Particle &pco, float (*op_fun)(const float *)) {
  for (uint32_t e = 0; e < pco.ev_count; e++) {
    update_up(pco, op_fun);
    std::cout << "EV: " << e << std::endl;
    if (pco.LOG_BEST_DATA)
      std::cout << "BEST: " << op_fun(pco.g_best_p) << std::endl;
    if (pco.LOG_AVERAGE_DATA)
      std::cout << "AVERAGE: "
                << average(pco.position_vector, op_fun, pco.po_size,
                           pco.pos_len)
                << std::endl;
    std::cout << "========================" << std::endl;
  }
}

float average(const float *position_vector, float (*op_fun)(const float *),
              uint32_t len1, uint32_t len2) {
  float sum = 0.0f;
  for (uint32_t i = 0; i < len1; i++) {
    sum += op_fun(position_vector + i * len2);
  }
  return sum / static_cast<float>(len1);
}

void __single_threads_down__(Particle &pco, float (*op_fun)(const float *),
                             uint32_t start_, uint32_t end_) {
  for (uint32_t p = start_; p < end_; p++) {
    for (uint32_t i = 0; i < pco.pos_len; i++) {
      pco.velocity_vector[p * pco.pos_len + i] =
          pco.inertia_weight * pco.velocity_vector[p * pco.pos_len + i] +
          pco.acc_factor1 * random(gen) *
              (pco.p_best_p[p * pco.pos_len + i] -
               pco.position_vector[p * pco.pos_len + i]) +
          pco.acc_factor2 * random(gen) *
              (pco.g_best_p[i] - pco.position_vector[p * pco.pos_len + i]);

      pco.position_vector[p * pco.pos_len + i] =
          pco.position_vector[p * pco.pos_len + i] +
          pco.velocity_vector[p * pco.pos_len + i];
    }
    pco.p_func_value[p] = op_fun(pco.position_vector + p * pco.pos_len);

    if (op_fun(pco.position_vector + p * pco.pos_len) <
        op_fun(pco.p_best_p + p * pco.pos_len)) {
      memcpy(pco.p_best_p + p * pco.pos_len,
             pco.position_vector + p * pco.pos_len,
             pco.pos_len * sizeof(float));
    }
  }
}

void update_down(Particle &pco, float (*op_fun)(const float *)) {

  uint32_t t_nums = (pco.po_size - 1ul) / pco.threads_t + 1;

  uint32_t start_;
  uint32_t end_;

  std::vector<std::thread> v_threads;
  v_threads.reserve(pco.threads_t);

  for (uint32_t t = 0; t < pco.threads_t; t++) {

    start_ = t * t_nums;
    end_ = (t + 1) * t_nums;
    if (end_ > pco.po_size)
      end_ = pco.po_size;
    v_threads.emplace_back(__single_threads_down__, pco, op_fun, start_, end_);
  }
  for (auto &th : v_threads)
    th.join();

  uint32_t g_best_index =
      find_best_position_down(pco.p_func_value, pco.po_size);
  memcpy(pco.g_best_p, pco.position_vector + g_best_index * pco.pos_len,
         pco.pos_len * sizeof(float));
}

void __single_threads_up__(Particle &pco, float (*op_fun)(const float *),
                           uint32_t start_, uint32_t end_) {
  for (uint32_t p = start_; p < end_; p++) {
    for (uint32_t i = 0; i < pco.pos_len; i++) {
      pco.velocity_vector[p * pco.pos_len + i] =
          pco.inertia_weight * pco.velocity_vector[p * pco.pos_len + i] +
          pco.acc_factor1 * random(gen) *
              (pco.p_best_p[p * pco.pos_len + i] -
               pco.position_vector[p * pco.pos_len + i]) +
          pco.acc_factor2 * random(gen) *
              (pco.g_best_p[i] - pco.position_vector[p * pco.pos_len + i]);

      pco.position_vector[p * pco.pos_len + i] =
          pco.position_vector[p * pco.pos_len + i] +
          pco.velocity_vector[p * pco.pos_len + i];
    }
    pco.p_func_value[p] = op_fun(pco.position_vector + p * pco.pos_len);
    if (op_fun(pco.position_vector + p * pco.pos_len) >
        op_fun(pco.p_best_p + p * pco.pos_len)) {
      memcpy(pco.p_best_p + p * pco.pos_len,
             pco.position_vector + p * pco.pos_len,
             pco.pos_len * sizeof(float));
    }
  }
}

void update_up(Particle &pco, float (*op_fun)(const float *)) {
  uint32_t t_nums = (pco.po_size - 1ul) / pco.threads_t + 1;

  uint32_t start_;
  uint32_t end_;

  std::vector<std::thread> v_threads;
  v_threads.reserve(pco.threads_t);

  for (uint32_t t = 0; t < pco.threads_t; t++) {
    start_ = t * t_nums;
    end_ = (t + 1) * t_nums;
    if (end_ > pco.po_size)
      end_ = pco.po_size;
    v_threads.emplace_back(__single_threads_up__, pco, op_fun, start_, end_);
  }
  for (auto &th : v_threads)
    th.join();
  uint32_t g_best_index = find_best_position_up(pco.p_func_value, pco.po_size);
  memcpy(pco.g_best_p, pco.position_vector + g_best_index * pco.pos_len,
         pco.pos_len * sizeof(float));
}

void destroy_particle(Particle &pco) {
  free(pco.velocity_vector);
  free(pco.position_vector);
  free(pco.p_best_p);
  free(pco.g_best_p);
  free(pco.p_func_value);
}

} // namespace particle
