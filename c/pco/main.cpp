// Copyright by 2022.3 3227068950@qq.com
// author 3227068950@qq.com

#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <stdint.h>

#define OPTIMIZE_UP 0

#define TEST 0

static uint32_t po_size = 100ul;
static uint32_t ev_count = 100ul;

static uint32_t pos_len = 12ul;

static float inertia_weight = 0.4f;
static float acc_factor1 = 2.f;
static float acc_factor2 = 2.f;

static std::random_device rd;
static std::mt19937 gen;
static std::uniform_real_distribution<float> random;

struct particle_res {
  uint32_t best_index = 0ul;
  float *best_res = nullptr;
};

float op_func(const float *arr);

void define_solution_domain(float *arr);

void define_velocity_range(float *arr);

void init_particle(float *velocity_vector, float *position_vector,
                   float *p_best_p, float *g_best_p);

uint32_t find_best_position(float (*op_fun)(const float *), const float *arr);

uint32_t start(float *velocity_vector, float *position_vector, float *p_best_p,
               float *g_best_p);

void update_vp(float *velocity_vector, float *position_vector,
               const float *p_best_p, const float *g_best_p);

void update_fitness(const float *position_vector, float *p_best_p,
                    float *g_best_p);

void show_array(const float *arr, uint32_t len1, uint32_t len2);

float func_average(const float *position_vector);

int main() {
  auto velocity_vector =
      static_cast<float *>(malloc(po_size * pos_len * sizeof(float)));
  auto position_vector =
      static_cast<float *>(malloc(po_size * pos_len * sizeof(float)));
  auto p_best_p =
      static_cast<float *>(malloc(po_size * pos_len * sizeof(float)));

  auto g_best_p = static_cast<float *>(malloc(pos_len * sizeof(float)));

  init_particle(velocity_vector, position_vector, p_best_p, g_best_p);
  
#if TEST
  auto func_arr = static_cast<float *>(malloc(po_size * sizeof(float)));

  for (uint32_t i = 0; i < po_size; i++) {
    func_arr[i] = op_func(position_vector + i * pos_len);
  }

  std::cout << "\nvelocity vector:" << std::endl;
  show_array(velocity_vector, po_size, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\nposition vector:" << std::endl;
  show_array(position_vector, po_size, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\npersonal best position: " << std::endl;
  show_array(p_best_p, po_size, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\nglobal best position: " << std::endl;
  show_array(g_best_p, 1, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\nfunc value: " << std::endl;
  show_array(func_arr, po_size, 1);
  free(func_arr);
#endif

  start(velocity_vector, position_vector, p_best_p, g_best_p);

  show_array(g_best_p, 1, pos_len);

  free(velocity_vector);
  free(position_vector);
  free(p_best_p);
  free(g_best_p);
}

float op_func(const float *arr) {
  float sum = 0.f;
  for (uint32_t i = 0; i < pos_len; i++) {
    sum += (arr[i] - static_cast<float>(i)) * (arr[i] - static_cast<float>(i));
  }
  return sum;
}

void init_particle(float *velocity_vector, float *position_vector,
                   float *p_best_p, float *g_best_p) {
  for (uint32_t i = 0; i < po_size; i++) {
    define_solution_domain(position_vector + i * pos_len);
    define_velocity_range(velocity_vector + i * pos_len);
  }
  memcpy(p_best_p, position_vector, po_size * pos_len * sizeof(float));
  uint32_t g_best_index = find_best_position(op_func, position_vector);
  memcpy(g_best_p, position_vector + g_best_index * pos_len,
         pos_len * sizeof(float));
}

void define_solution_domain(float *arr) {
  for (uint32_t i = 0; i < pos_len; i++) {
    arr[i] = 30 * (random(gen) - 0.5f);
  }
}

void define_velocity_range(float *arr) { define_solution_domain(arr); }

uint32_t find_best_position(float (*op_fun)(const float *), const float *arr) {
  uint32_t index = 0ul;
  float best_value = op_fun(arr);
  float temp_value;
  for (uint32_t i = 1; i < po_size; i++) {
    temp_value = op_fun(arr + i * pos_len);
#if OPTIMIZE_UP
    if (best_value < temp_value) {
      best_value = temp_value;
      index = i;
    }
#else
    if (best_value > temp_value) {
      best_value = temp_value;
      index = i;
    }
#endif
  }
  return index;
}

void show_array(const float *arr, uint32_t len1, uint32_t len2) {
  for (uint32_t i = 0; i < len1; i++) {
    for (uint32_t j = 0; j < len2; j++) {
      std::cout << arr[i * len2 + j] << " ";
    }
    std::cout << std::endl;
  }
}

uint32_t start(float *velocity_vector, float *position_vector, float *p_best_p,
               float *g_best_p) {

  for (uint32_t e = 0; e < ev_count; e++) {
    update_vp(velocity_vector, position_vector, p_best_p, g_best_p);
    update_fitness(position_vector, p_best_p, g_best_p);
    std::cout << "EV: " << e << std::endl;
    std::cout << "BEST: " << op_func(g_best_p) << std::endl;
    std::cout << "AVERAGE: " << func_average(position_vector) << std::endl;
    std::cout << "========================" << std::endl;
  }
  return 0ul;
}

void update_vp(float *velocity_vector, float *position_vector,
               const float *p_best_p, const float *g_best_p) {
  for (uint32_t p = 0; p < po_size; p++) {
    for (uint32_t i = 0; i < pos_len; i++) {
      velocity_vector[p * pos_len + i] =
          inertia_weight * velocity_vector[p * pos_len + i] +
          acc_factor1 * random(gen) *
              (p_best_p[p * pos_len + i] - position_vector[p * pos_len + i]) +
          acc_factor2 * random(gen) *
              (g_best_p[i] - position_vector[p * pos_len + i]);

      position_vector[p * pos_len + i] =
          position_vector[p * pos_len + i] + velocity_vector[p * pos_len + i];
    }
  }
}

void update_fitness(const float *position_vector, float *p_best_p,
                    float *g_best_p) {
  for (uint32_t i = 0; i < po_size; i++) {
#if OPTIMIZE_UP
    if (op_func(position_vector + i * pos_len) >
        op_func(p_best_p + i * pos_len)) {
      memcpy(p_best_p + i * pos_len, position_vector + i * pos_len,
             pos_len * sizeof(float));
    }
#else
    if (op_func(position_vector + i * pos_len) <
        op_func(p_best_p + i * pos_len)) {
      memcpy(p_best_p + i * pos_len, position_vector + i * pos_len,
             pos_len * sizeof(float));
    }
#endif
  }
  uint32_t g_best_index = find_best_position(op_func, position_vector);
  memcpy(g_best_p, position_vector + g_best_index * pos_len,
         pos_len * sizeof(float));
}
float func_average(const float *position_vector) {
  float sum = 0.0f;
  for (uint32_t i = 0; i < po_size; i++) {
    sum += op_func(position_vector + i * pos_len);
  }
  return sum / static_cast<float>(po_size);
}
