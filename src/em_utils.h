// ----------------------------------------------------------------------------
//
// Copyright (C) 2023 by The authors
//
// This file is part of the BSEM3DFWD program.
//
// This code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.
//
// ----------------------------------------------------------------------------

#ifndef _EM_UTILS_H_
#define _EM_UTILS_H_ 1

#include <petsc.h>

#include <string>

std::string string_format(const char *, ...);
std::string parse_string(const std::string &);

/**
 * @brief A helper class for logging events using PetscLogEvent.
 *
 * This class provides a convenient way to log events using PetscLogEvent.
 * It automatically calls PetscLogEventBegin() and PetscLogEventEnd()
 * in its constructor and destructor, respectively.
 */
class LogEventHelper {
public:
  LogEventHelper(PetscLogEvent ple) : ple_(ple) { PetscLogEventBegin(ple_, 0, 0, 0, 0); }
  ~LogEventHelper() { PetscLogEventEnd(ple_, 0, 0, 0, 0); }

private:
  PetscLogEvent ple_;
};

/**
 * @brief A helper class for logging stages using PetscLogStage.
 *
 * This class provides a convenient way to log stages using PetscLogStage.
 * It automatically calls PetscLogStageRegister(), PetscLogStagePush(),
 * and PetscLogStagePop() in its constructor and destructor, respectively.
 */
class LogStageHelper {
public:
  LogStageHelper(const std::string &name) {
    PetscLogStage pls;
    PetscLogStageRegister(name.c_str(), &pls);
    PetscLogStagePush(pls);
  }
  ~LogStageHelper() { PetscLogStagePop(); }
};

#endif

