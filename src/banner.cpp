#include "banner.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>

#ifdef _WIN32
  #include <io.h>
  #define isatty _isatty
  #define fileno _fileno
#else
  #include <unistd.h>
#endif

#include "generated/version.h"

namespace amdat { namespace cli {

static std::string detect_platform() {
#if defined(__linux__)
  return "x86_64-linux-gnu";
#elif defined(__APPLE__)
  return "apple-darwin";
#elif defined(_WIN32)
  return "win64";
#else
  return "unknown";
#endif
}

bool stdout_is_tty() { return isatty(fileno(stdout)); }

bool env_no_banner() {
  const char* v = std::getenv("AMDAT_NO_BANNER");
  return v && std::string(v) == "1";
}

VersionInfo collect_version_info() {
  VersionInfo v;
  v.name         = amdat::build::NAME;
  v.semver       = amdat::build::SEMVER;
  v.git_commit   = amdat::build::GIT_COMMIT;
  v.git_describe = amdat::build::GIT_DESCRIBE;
  v.build_date   = amdat::build::BUILD_DATE;
  v.compiler     = amdat::build::COMPILER;
  v.openmp       = amdat::build::OPENMP;
  v.platform     = detect_platform();
  return v;
}

void print_one_line_banner_stderr(const VersionInfo& b) {
  // Note: stderr keeps stdout clean for pipelines
  std::fprintf(stderr, "%s %s (commit %s, built %s) — Advanced Molecular Dynamics Analysis Toolkit\n",
               b.name.c_str(), b.semver.c_str(), b.git_commit.c_str(), b.build_date.c_str());
}

void print_one_line_banner(const VersionInfo& b) {
  std::cout << b.name << ' ' << b.semver
            << " (commit " << b.git_commit
            << ", built " << b.build_date
            << ") - Amorphous Molecular Dynamics Analysis Toolkit\n";
}


void print_version_text(const VersionInfo& b) {
  std::printf(
    "%s %s\n"
    "commit: %s\n"
    "describe: %s\n"
    "built: %s\n"
    "compiler: %s\n"
    "openmp: %s\n"
    "platform: %s\n"
    "license: GPL-3.0-or-later\n"
    "homepage: https://dssimmons-codes.github.io/AMDAT/\n"
    "cite: Please cite AMDAT (see docs for DOI)\n",
    b.name.c_str(), b.semver.c_str(),
    b.git_commit.c_str(),
    b.git_describe.c_str(),
    b.build_date.c_str(),
    b.compiler.c_str(),
    b.openmp ? "on" : "off",
    b.platform.c_str()
  );
}

void print_version_json(const VersionInfo& b) {
  std::printf(
    "{\n"
    "  \"name\": \"%s\",\n"
    "  \"version\": \"%s\",\n"
    "  \"git\": {\"commit\": \"%s\", \"describe\": \"%s\"},\n"
    "  \"build\": {\"date\": \"%s\", \"compiler\": \"%s\", \"openmp\": %s},\n"
    "  \"platform\": \"%s\",\n"
    "  \"license\": \"GPL-3.0-or-later\",\n"
    "  \"homepage\": \"https://dssimmons-codes.github.io/AMDAT/\"\n"
    "}\n",
    b.name.c_str(), b.semver.c_str(),
    b.git_commit.c_str(), b.git_describe.c_str(),
    b.build_date.c_str(), b.compiler.c_str(), (b.openmp ? "true" : "false"),
    b.platform.c_str()
  );
}

}} // namespace amdat::cli
