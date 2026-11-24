// Read-only detection of the CPU and memory this process may use (plink4r_resources() in R).
//
// Common setups only, by design:
// threads: the scheduler's allocation (SLURM_CPUS_PER_TASK, SLURM_CPUS_ON_NODE; PBS/Torque NCPUS, PBS_NP;
// LSF LSB_DJOB_NUMPROC; SGE NSLOTS), capped by the affinity mask; with no scheduler, the affinity mask when it
// is narrower than the machine (a cgroup / taskset limit), else 4. The total core count of a node is never
// used: on a shared login node or workstation it belongs to everyone. Windows / macOS: 4.
// mem_gb: half of the smallest of the SLURM allocation (SLURM_MEM_PER_NODE, or SLURM_MEM_PER_CPU x threads),
// the cgroup limit (v2 memory.max, v1 memory.limit_in_bytes) and MemAvailable; with neither a scheduler
// allocation nor a cgroup limit (a shared machine) at most 16 GB; 2 GB when nothing is known.
// Nothing here writes an environment variable or any global thread setting.
#ifndef PLINK4R_RESOURCES_H
#define PLINK4R_RESOURCES_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#if defined(__linux__)
#include <sched.h>
#include <unistd.h>
#endif

namespace plink4r {

struct Resources {
  int threads = 4;
  double mem_gb = 2;
  std::string threads_from = "default";
  std::string mem_from = "default";
};

namespace res_detail {

inline long env_long(const char* name) {
  const char* v = std::getenv(name);
  if (!v || !*v) return -1;
  char* end = nullptr;
  const long x = std::strtol(v, &end, 10);
  return (end != v && x > 0) ? x : -1;
}

// SLURM memory strings are in MB unless suffixed (K, M, G, T).
inline double env_mb(const char* name) {
  const char* v = std::getenv(name);
  if (!v || !*v) return -1;
  char* end = nullptr;
  const double x = std::strtod(v, &end);
  if (end == v || !(x > 0)) return -1;
  switch (*end) {
    case 'K': case 'k': return x / 1024;
    case 'G': case 'g': return x * 1024;
    case 'T': case 't': return x * 1024 * 1024;
    default: return x;
  }
}

inline double file_bytes(const std::string& path) {
  std::ifstream f(path);
  if (!f) return -1;
  std::string s;
  f >> s;
  if (s.empty() || s == "max") return -1;
  const double x = std::strtod(s.c_str(), nullptr);
  return (x > 0 && x < 1e18) ? x : -1;  // v1 reports ~2^63 for "no limit"
}

inline double cgroup_mb() {
#if defined(__linux__)
  std::ifstream f("/proc/self/cgroup");
  std::string line, v2path, v1path;
  while (std::getline(f, line)) {
    const auto a = line.find(':'), b = line.find(':', a + 1);
    if (a == std::string::npos || b == std::string::npos) continue;
    const std::string ctl = line.substr(a + 1, b - a - 1), path = line.substr(b + 1);
    if (line.compare(0, 2, "0:") == 0 && ctl.empty()) v2path = path;
    if (ctl.find("memory") != std::string::npos) v1path = path;
  }
  double best = -1;
  auto take = [&](double x) { if (x > 0 && (best < 0 || x < best)) best = x; };
  take(file_bytes("/sys/fs/cgroup" + v2path + "/memory.max"));
  take(file_bytes("/sys/fs/cgroup/memory.max"));
  take(file_bytes("/sys/fs/cgroup/memory" + v1path + "/memory.limit_in_bytes"));
  return best > 0 ? best / (1024.0 * 1024.0) : -1;
#else
  return -1;
#endif
}

inline double meminfo_mb() {
#if defined(__linux__)
  std::ifstream f("/proc/meminfo");
  std::string key;
  double kb;
  std::string unit;
  while (f >> key >> kb) {
    std::getline(f, unit);
    if (key == "MemAvailable:") return kb / 1024.0;
  }
#endif
  return -1;
}

inline int affinity_cpus(bool* narrower) {
  *narrower = false;
#if defined(__linux__) && defined(CPU_COUNT)
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) != 0) return -1;
  const int n = CPU_COUNT(&set);
  const long online = sysconf(_SC_NPROCESSORS_ONLN);
  *narrower = online > 0 && n < online;
  return n > 0 ? n : -1;
#else
  return -1;
#endif
}

}  // namespace res_detail

// threads <= 0 or mem_gb <= 0: detect; a positive value is used as given.
inline Resources detect_resources(int threads = 0, double mem_gb = 0) {
  using namespace res_detail;
  Resources r;
  if (threads > 0) {
    r.threads = threads; r.threads_from = "argument";
  } else {
    bool narrower = false;
    const int aff = affinity_cpus(&narrower);
    static const char* const sched[] = {"SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "NCPUS", "PBS_NP",
                                        "LSB_DJOB_NUMPROC", "NSLOTS"};
    long s = -1;
    std::string from;
    for (const char* v : sched) { s = env_long(v); if (s > 0) { from = v; break; } }
    if (s > 0) {
      r.threads = static_cast<int>(aff > 0 ? std::min<long>(s, aff) : s); r.threads_from = from;
    } else if (aff > 0 && narrower) {
      r.threads = aff; r.threads_from = "affinity";
    }
  }
  if (mem_gb > 0) {
    r.mem_gb = mem_gb; r.mem_from = "argument";
  } else {
    double best = -1;
    auto take = [&](double mb, const char* from) { if (mb > 0 && (best < 0 || mb < best)) { best = mb; r.mem_from = from; } };
    const double node = env_mb("SLURM_MEM_PER_NODE");
    if (node > 0) take(node, "SLURM_MEM_PER_NODE");
    else { const double cpu = env_mb("SLURM_MEM_PER_CPU"); if (cpu > 0) take(cpu * r.threads, "SLURM_MEM_PER_CPU"); }
    const double cg = cgroup_mb();
    take(cg, "cgroup");
    const bool limited = best > 0;  // a scheduler allocation or a cgroup limit
    take(meminfo_mb(), "MemAvailable");
    if (best > 0) r.mem_gb = best / 2048.0;  // half, in GB
    if (!limited && r.mem_gb > 16) { r.mem_gb = 16; r.mem_from = "MemAvailable (capped at 16 GB on a shared machine)"; }
  }
  return r;
}

inline int resolve_threads(int threads) { return threads > 0 ? threads : detect_resources(threads, 1).threads; }
inline double resolve_mem_gb(double mem_gb, int threads) { return mem_gb > 0 ? mem_gb : detect_resources(threads, mem_gb).mem_gb; }

}  // namespace plink4r
#endif
