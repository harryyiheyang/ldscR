// Run-time choice of the SIMD matrix-product kernel (la_kernel.cpp, compiled for generic, AVX2+FMA and AVX-512
// by src/Makevars), as in plink4r: R's package flags never define __AVX2__, so Eigen products in ordinary
// translation units run on SSE2. la() returns the best table the CPU supports, chosen once per process;
// LDSCR_SIMD=avx512|avx2|scalar caps the choice.
#ifndef LDSCR_LA_DISPATCH_H
#define LDSCR_LA_DISPATCH_H

#include <cstdlib>
#include <string>

namespace ldscr {

struct LaFns {
  // Out (M x N, ldo) = A (M x K, lda) * B (K x N, ldb); column-major floats, single thread.
  void (*sgemm_nn)(const float* A, int M, int K, int lda, const float* B, int N, int ldb, float* Out, int ldo);
  const char* isa;
};

}  // namespace ldscr

extern const ldscr::LaFns ldscr_la_generic;
#if defined(LDSCR_LA_X86)
extern const ldscr::LaFns ldscr_la_avx2;
extern const ldscr::LaFns ldscr_la_avx512;
#endif

namespace ldscr {

inline const LaFns& pick_la() {
  const char* env = std::getenv("LDSCR_SIMD");
  const std::string want = env ? env : "";
#if defined(LDSCR_LA_X86) && (defined(__GNUC__) || defined(__clang__))
  __builtin_cpu_init();
  const bool a512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
                    __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw");
  const bool a2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
  if (a512 && (want.empty() || want == "avx512")) return ldscr_la_avx512;
  if (a2 && (want.empty() || want == "avx512" || want == "avx2")) return ldscr_la_avx2;
#endif
  return ldscr_la_generic;
}
inline const LaFns& la() { static const LaFns& t = pick_la(); return t; }

}  // namespace ldscr
#endif
