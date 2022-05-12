

#include <iostream>
#include <ctime>

#define OPTIMIZE_UP 0

#include "particle_v2.h"

#define TEST 0

static uint32_t po_size = 100000ul;
static uint32_t ev_count = 100ul;

static uint32_t pos_len = 13ul;

static float inertia_weight = 0.3f;
static float acc_factor1 = 2.f;
static float acc_factor2 = 2.f;

static uint32_t threads_num = 20;

float op_func(const float *arr);

void define_solution_domain(float *arr);

void define_velocity_range(float *arr);

int main() {
  particle::Particle pco;
  pco.po_size = po_size;
  pco.pos_len = pos_len;
  pco.ev_count = ev_count;
  pco.inertia_weight = inertia_weight;
  pco.acc_factor1 = acc_factor1;
  pco.acc_factor2 = acc_factor2;
  pco.threads_t = threads_num;

  // pco.LOG_AVERAGE_DATA = false;
  pco.thread_api = particle::STD_THREAD;

  particle::init_particle(pco, op_func, define_solution_domain,
                          define_velocity_range);
#if TEST
  auto func_arr = static_cast<float *>(malloc(po_size * sizeof(float)));

  for (uint32_t i = 0; i < po_size; i++) {
    func_arr[i] = op_func(pco.position_vector + i * pos_len);
  }

  std::cout << "\nvelocity vector:" << std::endl;
  particle::show_array(pco.velocity_vector, po_size, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\nposition vector:" << std::endl;
  particle::show_array(pco.position_vector, po_size, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\npersonal best position: " << std::endl;
  particle::show_array(pco.p_best_p, po_size, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\nglobal best position: " << std::endl;
  particle::show_array(pco.g_best_p, 1, pos_len);
  std::cout << "==========================" << std::endl;
  std::cout << "\nfunc value: " << std::endl;
  particle::show_array(func_arr, po_size, 1);
  free(func_arr);
#endif

  uint32_t start_time = std::clock();
  particle::start(pco, op_func);
  uint32_t stop_time = std::clock();
  
  particle::show_array(pco.g_best_p, 1, pos_len);

  std::cout << "Elapsed time: " << stop_time - start_time << "ms" << " for threads " << pco.threads_t << std::endl;
  particle::destroy_particle(pco);
}

float op_func(const float *arr) {
  float sum = 0.f;
  for (uint32_t i = 0; i < pos_len; i++) {
    sum += (arr[i] - static_cast<float>(i)) * (arr[i] - static_cast<float>(i));
  }
  return sum + 1.f;
}

void define_solution_domain(float *arr) {
  for (uint32_t i = 0; i < pos_len; i++) {
    arr[i] = 30 * (particle::random(particle::gen) - 0.5f);
  }
}

void define_velocity_range(float *arr) { define_solution_domain(arr); }
