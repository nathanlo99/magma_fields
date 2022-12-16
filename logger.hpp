
#pragma once

#include "timing.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

enum class LogFile { None, Stderr };

inline std::ostream &log(const LogFile log = LogFile::Stderr) {
  static std::ostringstream dev_null;
  switch (log) {
  case LogFile::Stderr:
    return std::cerr << "LOG [" << std::setw(6) << get_elapsed_ms() << "ms]: ";
  default:
    return dev_null;
  }
  __builtin_unreachable();
}
