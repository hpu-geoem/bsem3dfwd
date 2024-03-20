CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(bsem3dfwd CXX)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

FIND_PACKAGE(deal.II 9.5.0 QUIET HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a sufficiently recent version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT DEAL_II_WITH_PETSC)
  MESSAGE(FATAL_ERROR "\n"
    "*** deal.II MUST be configured with PETSc. ***\n"
    )
ENDIF()

IF (${DEAL_II_PETSC_VERSION} VERSION_LESS 3.18.0)
  MESSAGE(FATAL_ERROR "\n"
    "*** deal.II MUST be configured with PETSc 3.18.0 or higher. ***\n"
    )
ENDIF()

IF (NOT DEAL_II_PETSC_WITH_HYPRE)
  MESSAGE(FATAL_ERROR "\n"
    "*** deal.II MUST be configured with PETSc with HYPRE to use the AMS preconditioner. ***\n"
    )
ENDIF()

IF (DEAL_II_PETSC_WITH_COMPLEX)
  MESSAGE(FATAL_ERROR "\n"
    "*** deal.II MUST be configured with PETSc not with complex scalars. ***\n"
    )
ENDIF()

SET(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER})
SET(CMAKE_CXX_COMPILER_FLAGS ${DEAL_II_CXX_COMPILER_FLAGS})

INCLUDE_DIRECTORIES(${DEAL_II_INCLUDE_DIRS})

ADD_SUBDIRECTORY(src)
