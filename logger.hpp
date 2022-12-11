
#pragma once

#include <fstream>
#include <iostream>
#include <sstream>

enum LogFile { None, Stderr };

inline std::ostream &log(const LogFile log = Stderr) {
  static std::ostringstream dev_null;
  switch (log) {
  case Stderr:
    return std::cerr;
  default:
    return dev_null;
  }
  __builtin_unreachable();
}
