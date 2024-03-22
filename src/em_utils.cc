// ----------------------------------------------------------------------------
//
// Copyright (C) 2023 by Ce Qin <ce.qin@hpu.edu.cn>
//
// This file is part of the BSEM3DFWD program.
//
// This code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.
//
// ----------------------------------------------------------------------------

#include <stdarg.h>
#include <string.h>

#include <cmath>
#include <memory>
#include <string>

/**
 * @brief Formats a string using a printf-like syntax.
 *
 * This function formats a string using a printf-like syntax and returns the
 * resulting string. The function takes a format string followed by a variable
 * number of arguments, which are used to replace the format specifiers in the
 * format string.
 *
 * @param fmt_str The format string.
 * @param ... The variable arguments to replace the format specifiers.
 *
 * @return The formatted string.
 */
std::string string_format(const char *fmt_str, ...) {
  va_list ap;
  int final_n, n;
  std::string str;
  std::unique_ptr<char[]> formatted;

  n = strlen(fmt_str) * 2;

  while (1) {
    formatted.reset(new char[n]);
    strcpy(&formatted[0], fmt_str);

    va_start(ap, fmt_str);
    final_n = std::vsnprintf(&formatted[0], n, fmt_str, ap);
    va_end(ap);

    if (final_n < 0 || final_n >= n) {
      n += std::abs(final_n - n + 1);
    } else {
      break;
    }
  }

  return std::string(formatted.get());
}

/**
 * @brief Parses a string by removing any trailing or leading whitespace and
 *        any characters after a '#' character.
 *
 * This function takes a string and removes any trailing or leading whitespace
 * characters. It also removes any characters after the first occurrence of a
 * '#' character. The resulting string is returned.
 *
 * @param s The input string to parse.
 *
 * @return The parsed string.
 */
std::string parse_string(const std::string &s) {
  size_t pos;
  std::string l = s;

  pos = l.find("#");
  if (pos != std::string::npos) {
    l.erase(pos, std::string::npos);
  }

  pos = l.find_first_not_of(" \t");
  if (pos != 0 && pos != std::string::npos) {
    l.erase(l.begin(), l.begin() + pos);
  }

  pos = l.find_last_not_of(" \t");
  if (pos != (l.size() - 1) && pos != std::string::npos) {
    l.erase(pos + 1, std::string::npos);
  }

  if (l.size() > 0) {
    l += "\n";
  }

  return l;
}
