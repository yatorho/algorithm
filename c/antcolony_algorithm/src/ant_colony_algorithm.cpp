// Copyright by 2022.3 3227068950@qq.com
// author yatorho

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <random>
#include <stdint.h>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <windows.h>

#define NOT_IMPLEMENT_ERROR                                                    \
  std::cout << "NOT IMPLEMENT ERROR!" << std::endl;                            \
  std::exit(EXIT_FAILURE)

#define TEST 0
#define TEST_DATA_OUTPUT 0

#define USE_ANT_CYCLE_SYSTEM 1
#define USE_ANT_QUANTITY_SYSTEM 0
#define USE_ANT_DENSITY_SYSTEM 0

#define CHECK 0

#define THREAD_WIN 1

static uint32_t ant_count = 80;
static float heuristic_factor = 2.3f;
static float pheromone_factor = 1.2f;
static float pher_reduce_factor = 0.1f;
static const float pheromone_sum = 30.f;
static const uint32_t iteration = 100;
static float pher_default_value = 10.f;
static uint32_t city_count = 5;

static uint32_t thread_nums = 1;

void init_city_test(float *arr, uint32_t nums);

void show_array(const float *arr, uint32_t axis1, uint32_t axis2);

void init(const float *city, float *heur, float *pher);

uint32_t start_single(std::vector<uint32_t> *record, const float *city,
                      float *heur, float *pher);

uint32_t start_multi(std::vector<uint32_t> *record, const float *city,
                     float *heur, float *pher, uint32_t t_nums);

void __single_thread_f__(std::vector<uint32_t> *record, const float *heur,
                         float *pher, float *prob, const float *city,
                         float *path_len, const uint32_t &nums,
                         uint32_t thread_id, std::mt19937 &gen,
                         std::uniform_int_distribution<uint32_t> &city_dis,
                         std::uniform_real_distribution<float> &ran_dis);

DWORD WINAPI __single_thread_f2__(PVOID parm);

struct thread_parm {
  std::vector<uint32_t> *record = nullptr;
  const float *heur = nullptr;
  float *pher = nullptr;
  float *prob = nullptr;
  const float *city = nullptr;
  float *path_len = nullptr;
  uint32_t nums = 0;
  uint32_t thread_id = 0;
  std::mt19937 *gen = nullptr;
  std::uniform_int_distribution<uint32_t> *city_dis = nullptr;
  std::uniform_real_distribution<float> *ran_dis = nullptr;
};

uint32_t __roulette_select(const std::vector<uint32_t> &record,
                           std::uniform_real_distribution<float> &dis,
                           std::mt19937 &gen, const float *heur,
                           const float *pher, const float *prob);

void update_prob(const std::vector<uint32_t> &record, const float *heur,
                 const float *pher, float *prob);

void update_pher(const std::vector<uint32_t> *record, float *pher,
                 const float *path_len);

float asum_path_length(const std::vector<uint32_t> &record, const float *city);

float average(const float *arr, uint32_t len);

uint32_t index_min_array(const float *arr, uint32_t len);

int main() {
#if TEST
  auto city =
      static_cast<float *>(malloc(city_count * city_count * sizeof(float)));
  auto heur =
      static_cast<float *>(malloc(city_count * city_count * sizeof(float)));
  auto pher =
      static_cast<float *>(malloc(city_count * city_count * sizeof(float)));

  auto prob =
      static_cast<float *>(malloc(ant_count * city_count * sizeof(float)));

  auto record = new std::vector<uint32_t>[ant_count];

  init_city_test(city, city_count);
  init(city, heur, pher);
  show_array(city, city_count, city_count);
  show_array(pher, city_count, city_count);
  show_array(heur, city_count, city_count);
  std::cout << "INIT PASS!" << std::endl;

  record[0].push_back(2);
  record[0].push_back(1);

  record[2].push_back(4);
  record[6].push_back(3);
  record[6].push_back(0);

  update_prob(record[0], heur, pher, prob);
  update_prob(record[2], heur, pher, prob + 2 * city_count);
  update_prob(record[6], heur, pher, prob + 6 * city_count);
  show_array(prob, ant_count, city_count);
  std::cout << "UPDATE PROB PASS!" << std::endl;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis(0, 1);

  uint32_t res1 = __roulette_select(record[0], dis, gen, heur, pher, prob);
  uint32_t res2 =
      __roulette_select(record[6], dis, gen, heur, pher, prob + 6 * city_count);
  std::cout << "selected1: " << res1 << std::endl;
  std::cout << "selected2: " << res2 << std::endl;

  std::cout << "ROULETTE SELECT PASS!" << std::endl;

  uint32_t best_index;
  best_index = start_single(record, city, heur, pher);
  std::cout << "PATH No: " << best_index << ">>> ";

  for (uint32_t n = 0; n < city_count; n++) {
    if (n != city_count - 1)
      std::cout << record[best_index][n] << "->";
    else
      std::cout << record[best_index][n];
  }
  std::cout << std::endl;

  free(city);
  free(heur);
  free(pher);
  delete[] record;
  std::cout << "ALL PASS!" << std::endl;

#endif
  auto city =
      static_cast<float *>(malloc(city_count * city_count * sizeof(float)));
  auto heur =
      static_cast<float *>(malloc(city_count * city_count * sizeof(float)));
  auto pher =
      static_cast<float *>(malloc(city_count * city_count * sizeof(float)));

  auto prob =
      static_cast<float *>(malloc(ant_count * city_count * sizeof(float)));

  auto record = new std::vector<uint32_t>[ant_count];

  init_city_test(city, city_count);
  init(city, heur, pher);

  uint32_t best_index;
  // best_index = start_single(record, city, heur, pher);
  best_index = start_multi(record, city, heur, pher, thread_nums);
  std::cout << "PATH No: " << best_index << ">>> ";

  for (uint32_t n = 0; n < city_count; n++) {
    std::cout << record[best_index][n] << "->";
  }
  std::cout << record[best_index][0] << std::endl;

  free(city);
  free(heur);
  free(pher);
  delete[] record;
  std::cout << "ALGORITHM END!" << std::endl;
}

void init_city_test(float *arr, uint32_t nums) {
  if (nums != 5)
    std::cout
        << "WARNING! use other implementation to initialize city grah instead!"
        << std::endl;
  arr[0] = 1e-5f;
  arr[1] = 40.f;
  arr[2] = 55.f;
  arr[3] = 8.f;
  arr[4] = 30.f;
  arr[6] = 1e-5f;
  arr[7] = 15.f;
  arr[8] = 24.f;
  arr[9] = 62.f;
  arr[12] = 1e-5f;
  arr[13] = 20.f;
  arr[14] = 35.f;
  arr[18] = 1e-5f;
  arr[19] = 16.f;
  arr[24] = 1e-5f;
  for (uint32_t i = 0; i < nums; i++) {
    for (uint32_t j = 0; j < nums; j++) {
      if (i > j) {
        arr[i * nums + j] = arr[j * nums + i];
      }
    }
  }
}

void show_array(const float *arr, uint32_t axis1, uint32_t axis2) {
  for (uint32_t i = 0; i < axis1; i++) {
    for (uint32_t j = 0; j < axis2; j++) {
      std::cout << arr[i * axis2 + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void init(const float *city, float *heur, float *pher) {
  for (uint32_t i = 0; i < city_count * city_count; i++) {
    pher[i] = pher_default_value;
    heur[i] = 1.f / city[i];
  }
}

uint32_t start_single(std::vector<uint32_t> *record, const float *city,
                      float *heur, float *pher) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint32_t> city_dis(0, city_count - 1);
  std::uniform_real_distribution<float> ran_dis(0, 1);

  auto prob =
      static_cast<float *>(malloc(ant_count * city_count * sizeof(float)));

  auto path_len = static_cast<float *>(malloc(ant_count * sizeof(float)));

  for (uint32_t i = 0; i < iteration; i++) {

    for (uint32_t m = 0; m < ant_count; m++) {
      record[m].clear();
#if CHECK
      if (!record[m].empty()) {
        std::cout << "FATAL ERROR IN RECORD!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
#endif
      uint32_t start_city = city_dis(gen);
      record[m].push_back(start_city);

      uint32_t select;
      for (uint32_t c = 0; c < city_count - 1; c++) {
        update_prob(record[m], heur, pher, prob + m * city_count);
        select = __roulette_select(record[m], ran_dis, gen, heur, pher,
                                   prob + m * city_count);
#if CHECK
        if (select < 0 || select >= city_count) {
          std::cout << "FATAL ERROR IN __ROUTLETTE_SELECT!" << std::endl;
          std::cout << "select: " << select << std::endl;
          std::cout << "iteration: " << i << std::endl;
          std::cout << "antNo: " << m << std::endl;
          std::exit(EXIT_FAILURE);
        }
#endif
        record[m].push_back(select);
      }
#if CHECK
      if (record[m].size() != city_count) {
        std::cout << "FATAL ERROR IN RECORD ARRAY UPDATE!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
#endif
      path_len[m] = asum_path_length(record[m], city);
    }
    update_pher(record, pher, path_len);
    std::cout << "ITERATION: " << i << std::endl;
    std::cout << "BEST: " << path_len[index_min_array(path_len, ant_count)]
              << std::endl;
    std::cout << "AVERAGE: " << average(path_len, ant_count) << std::endl;
    std::cout << "====================================" << std::endl;
  }

  uint32_t best_index = index_min_array(path_len, ant_count);
  free(prob);
  free(path_len);

  return best_index; // output best path
}

uint32_t __roulette_select(const std::vector<uint32_t> &record,
                           std::uniform_real_distribution<float> &dis,
                           std::mt19937 &gen, const float *heur,
                           const float *pher, const float *prob) {
  float r = dis(gen);
  // float r = 0.5;
  float s = 0.0f;
#if TEST_DATA_OUTPUT
  std::cout << "random number: " << r << std::endl;
#endif
  for (uint32_t i = 0; i < city_count; i++) {
    s += prob[i];
#if TEST_DATA_OUTPUT
    std::cout << s << ", ";
#endif
    if (r < s) {
#if TEST_DATA_OUTPUT
      std::cout << std::endl;
#endif
#if CHECK
      if (i < 0 || i >= city_count) {
        std::cout << "FATAL ERROR IN __ROUTLETTE_SELECT!" << std::endl;
        std::cout << "selected: " << i << std::endl;
        std::exit(EXIT_FAILURE);
      }
#endif
      return i;
    }
  }

  std::cout << "FATAL ERROR>>> Failed to roulette selection!" << std::endl;
  std::exit(EXIT_FAILURE);
  return -1;
}

void update_prob(const std::vector<uint32_t> &record, const float *heur,
                 const float *pher, float *prob) {
  std::memset(prob, 0, city_count * sizeof(float));

  float prob_sum = 0.f;
  for (uint32_t c = 0; c < city_count; c++) {
    if (std::find(record.begin(), record.end(), c) == record.end()) {
      prob_sum += pow(pher[record.back() * city_count + c], pheromone_factor) *
                  pow(heur[record.back() * city_count + c], heuristic_factor);
#if CHECK
      if (std::isnan(prob_sum)) {
        std::cout << "WARNING>>> NAN!" << std::endl;
        // std::cout << "c: " << c << std::endl;
        // std::cout << "prob_sum: " << prob_sum << std::endl;
        std::cout << pher[record.back() * city_count + c] << std::endl;
        std::cout << heur[record.back() * city_count + c] << std::endl;
        std::cout << record.back() << std::endl;
        std::cout << record.back() << std::endl;
        std::cout << record.size() << std::endl;
      }
#endif
    }
  }
  for (uint32_t c = 0; c < city_count; c++) {
    if (std::find(record.begin(), record.end(), c) == record.end()) {
      prob[c] = pow(pher[record.back() * city_count + c], pheromone_factor) *
                pow(heur[record.back() * city_count + c], heuristic_factor) /
                prob_sum;
#if CHECK
      if (prob_sum == 0 || std::isnan(prob[c])) {
        std::cout << "FATAL ERROR IN UPDATE_PROB DIVIDED BY ZERO!" << std::endl;
        std::cout << "c: " << c << std::endl;
        std::cout << "prob_sum: " << prob_sum << std::endl;
        // std::cout << =
        std::exit(EXIT_FAILURE);
      }
#endif
    }
  }
}

void update_pher(const std::vector<uint32_t> *record, float *pher,
                 const float *path_len) {
  float delta_pher;
  for (uint32_t i = 0; i < city_count - 1; i++) {
    for (uint32_t j = i + 1; j < city_count; j++) {
      delta_pher = 0.f;

      for (uint32_t n = 0; n < ant_count; n++) {
        for (uint32_t c = 0; c < city_count - 1; c++) {
          if ((record[n][c] == i && record[n][c + 1] == j) ||
              (record[n][city_count - 1] == i && record[n][0] == j)) {
#if USE_ANT_CYCLE_SYSTEM
            delta_pher += pheromone_sum / path_len[n];
#elif USE_ANT_QUANTITY_SYSTEM
            NOT_IMPLEMENT_ERROR;
#else
            NOT_IMPLEMENT_ERROR;
#endif
          }
        }
      }
      pher[i * city_count + j] =
          (1 - pher_reduce_factor) * pher[i * city_count + j] + delta_pher;
    }
  }
  for (uint32_t i = 0; i < city_count; i++) {
    for (uint32_t j = 0; j <= i; j++) {
      pher[i * city_count + j] = pher[j * city_count + i];
    }
  }
#if CHECK
  for (uint32_t i = 0; i < city_count; i++) {
    for (uint32_t j = 0; j < city_count; j++) {
      if (pher[i * city_count + j] != pher[j * city_count + i]) {
        std::cout << "FATAL ERROR IN UPDATE PHREOMONE MATRIX" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
#endif
}

float asum_path_length(const std::vector<uint32_t> &record, const float *city) {
  float len = 0.f;
  for (uint32_t c = 0; c < city_count - 1; c++) {
#if CHECK
    if (record[c] < 0 || record[c] >= city_count || record[c + 1] < 0 ||
        record[c + 1] >= city_count) {
      std::cout << "FATAL ERROR IN ASUM PATH LENGTH!" << std::endl;
      std::cout << "record[" << c << "]: " << record[c] << std::endl;
      std::cout << "record[" << c + 1 << "]: " << record[c + 1] << std::endl;
      std::exit(EXIT_FAILURE);
    }
#endif
    len += city[record[c] * city_count + record[c + 1]];
  }
  len += city[record[city_count - 1] * city_count + record[0]];
  return len;
}

float average(const float *arr, uint32_t len) {
  float avg = 0.f;
  for (uint32_t i = 0; i < len; i++) {
    avg += arr[i];
  }
  return avg / len;
}

uint32_t index_min_array(const float *arr, uint32_t len) {
  uint32_t index = 0ul;
  for (uint32_t m = 0; m < ant_count; m++) {
    if (arr[index] > arr[m])
      index = m;
  }
  return index;
}

uint32_t start_multi(std::vector<uint32_t> *record, const float *city,
                     float *heur, float *pher, uint32_t t_nums) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint32_t> city_dis(0, city_count - 1);
  std::uniform_real_distribution<float> ran_dis(0, 1);

  const uint32_t nums = (ant_count - 1) / t_nums + 1;

  auto prob =
      static_cast<float *>(malloc(ant_count * city_count * sizeof(float)));

  auto path_len = static_cast<float *>(malloc(ant_count * sizeof(float)));

  for (uint32_t i = 0; i < iteration; i++) {
#if THREAD_WIN
    auto threads_parm = new thread_parm[t_nums];
    auto threads = static_cast<HANDLE *>(malloc(t_nums * sizeof(HANDLE)));
    for (uint32_t t = 0; t < t_nums; t++) {
      threads_parm[i].nums = nums;
      threads_parm[i].city = city;
      threads_parm[i].heur = heur;
      threads_parm[i].pher = pher;
      threads_parm[i].path_len = path_len;
      threads_parm[i].prob = prob;
      threads_parm[i].record = record;
      threads_parm[i].gen = &gen;
      threads_parm[i].city_dis = &city_dis;
      threads_parm[i].ran_dis = &ran_dis;
      threads_parm[i].thread_id = t;

      threads[i] = CreateThread(nullptr, 0, &__single_thread_f2__,
                                &threads_parm[i], 0, nullptr);
    }
    // std::cout << "success" << std::endl;
    WaitForMultipleObjects(t_nums, threads, TRUE, INFINITE);
#else
    {
      std::vector<std::thread> v_thread;
      v_thread.reserve(t_nums);
      for (uint32_t t = 0; t < t_nums; t++) {
        v_thread.emplace_back(&__single_thread_f__, record, heur, pher, prob,
                              city, path_len, nums, t, gen, city_dis, ran_dis);
      }
      for (auto &th : v_thread)
        th.join();
    }
#endif

    update_pher(record, pher, path_len);

    std::cout << "ITERATION: " << i << std::endl;
    std::cout << "BEST: " << path_len[index_min_array(path_len, ant_count)]
              << std::endl;
    std::cout << "AVERAGE: " << average(path_len, ant_count) << std::endl;
    std::cout << "====================================" << std::endl;
#if THREAD_WIN
    // delete[] threads_parm;
    // free(threads);
#endif
  }
  uint32_t best_index =
      index_min_array(path_len, ant_count); // output best path

  free(prob);
  free(path_len);

  return best_index; // output best path index
}

void __single_thread_f__(std::vector<uint32_t> *record, const float *heur,
                         float *pher, float *prob, const float *city,
                         float *path_len, const uint32_t &nums,
                         uint32_t thread_id, std::mt19937 &gen,
                         std::uniform_int_distribution<uint32_t> &city_dis,
                         std::uniform_real_distribution<float> &ran_dis) {
  uint32_t begin_ = nums * thread_id;
  uint32_t end_ = nums * (thread_id + 1);
  // uint32_t start;
  // uint32_t stop;

  if (end_ > ant_count)
    end_ = ant_count;

  // std::cout << "begin:" << begin_ << std::endl;
  // std::cout << "end: " << end_ << std::endl;
  // std::random_device rd;
  // std::mt19937 gen_(rd());
  // std::uniform_int_distribution<uint32_t> city_dis_(0, city_count - 1);
  // std::uniform_real_distribution<float> ran_dis_(0, 1);

  for (uint32_t n = begin_; n < end_; n++) {
    // if (n == (begin_ + 400))
    // start = std::clock();
    record[n].clear();

    uint32_t start_city = city_dis(gen);
    record[n].push_back(start_city);

    uint32_t select;

    for (uint32_t c = 0; c < city_count - 1; c++) {
      update_prob(record[n], heur, pher, prob + n * city_count);
      select = __roulette_select(record[n], ran_dis, gen, heur, pher,
                                 prob + n * city_count);
      record[n].push_back(select);
    }
    path_len[n] = asum_path_length(record[n], city);
  }
}

DWORD WINAPI __single_thread_f2__(PVOID parm) {
  auto p = static_cast<thread_parm *>(parm);

  uint32_t begin_ = p->nums * p->thread_id;
  uint32_t end_ = p->nums * (p->thread_id + 1);
  // std::cout << "begin:" << begin_ << std::endl;
  // std::cout << "end: " << end_ << std::endl;
  std::random_device rd;
  std::mt19937 gen_(rd());
  std::uniform_int_distribution<uint32_t> city_dis_(0, city_count - 1);
  std::uniform_real_distribution<float> ran_dis_(0, 1);

  if (end_ > ant_count)
    end_ = ant_count;

  for (uint32_t n = begin_; n < end_; n++) {
    p->record[n].clear();
    // uint32_t start_city = (*(p->city_dis))(*(p->gen));
    uint32_t start_city = city_dis_(gen_);
    p->record[n].push_back(start_city);

    uint32_t select;

    for (uint32_t c = 0; c < city_count - 1; c++) {
      update_prob(p->record[n], p->heur, p->pher, p->prob + n * city_count);
      select = __roulette_select(p->record[n], *(p->ran_dis), *(p->gen),
                                 p->heur, p->pher, p->prob + n * city_count);
      // select = __roulette_select(p->record[n], ran_dis_, gen_,
      //                            p->heur, p->pher, p->prob + n * city_count);
      //  std::cout << "SUCESS";
      p->record[n].push_back(select);
    }
    p->path_len[n] = asum_path_length(p->record[n], p->city);
  }

  return 0;
}
