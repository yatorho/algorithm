
#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <random>
#include <stdint.h>

#if OPTIMIZE_UP
#define start start_up
#define init_particle init_particle_up
#else
#define start start_down
#define init_particle init_particle_down
#endif

namespace particle {

enum THREAD_API { STD_THREAD = 0, WIN_THREAD = 1 };
struct Particle {
  uint32_t po_size = 100ul;
  uint32_t ev_count = 100ul;
  uint32_t pos_len = 12ul;

  float inertia_weight = 0.3f;
  float acc_factor1 = 2.f;
  float acc_factor2 = 2.f;

  float *velocity_vector = nullptr;
  float *position_vector = nullptr;
  float *p_best_p = nullptr;
  float *g_best_p = nullptr;

  float *p_func_value = nullptr;

  uint32_t threads_t = 1ul;

  bool LOG_AVERAGE_DATA = true;
  bool LOG_BEST_DATA = true;

  THREAD_API thread_api = STD_THREAD;
};

static std::mt19937 gen;
static std::uniform_real_distribution<float> random;

void init_particle_down(Particle &pco, float (*op_fun)(const float *),
                        void (*define_solution_domain)(float *),
                        void (*define_velocity_range)(float *));

void init_particle_up(Particle &pco, float (*op_fun)(const float *),
                      void (*define_solution_domain)(float *),
                      void (*define_velocity_range)(float *));

void destroy_particle(Particle &pco);

uint32_t find_best_position_down(float (*op_fun)(const float *),
                                 const float *arr, uint32_t po_size,
                                 uint32_t pos_len);
uint32_t find_best_position_up(float (*op_fun)(const float *), const float *arr,
                               uint32_t po_size, uint32_t pos_len);

void show_array(const float *arr, uint32_t len1, uint32_t len2);

void start_down(Particle &pco, float (*op_fun)(const float *));

void start_up(Particle &pco, float (*op_fun)(const float *));

float average(const float *position_vector, float (*op_fun)(const float *),
              uint32_t len1, uint32_t len2);

void update_down(Particle &pco, float (*op_fun)(const float *));
void update_up(Particle &pco, float (*op_fun)(const float *));

} // namespace particle
#endif // PARTICLE_H_
