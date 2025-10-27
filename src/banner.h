#pragma once
#include <string>

namespace amdat { namespace cli {

struct VersionInfo {
  std::string name;
  std::string semver;
  std::string git_commit;
  std::string git_describe;
  std::string build_date;
  std::string compiler;
  bool        openmp;
  std::string platform;
};

bool stdout_is_tty();
bool env_no_banner();
void print_one_line_banner_stderr(const VersionInfo& info);
void print_one_line_banner(const VersionInfo& info);
void print_version_text(const VersionInfo& info);
void print_version_json(const VersionInfo& info);
VersionInfo collect_version_info();

}} // namespace amdat::cli
