cmake_minimum_required(VERSION 3.20)
project(Forget21)

INCLUDE(CheckCXXSourceCompiles)

set(CPLEX_STUDIO
  "$ENV{CPLEX_1210}"
  CACHE STRING "The path to the top level CPLEX studio install directory"
)

include_directories(Code)
include_directories(${CPLEX_STUDIO}/cplex/include)
include_directories(${CPLEX_STUDIO}/concert/include)

add_compile_definitions(IL_STD)

set(source_files
  Code/BranchAndBound.cpp
  Code/Hyperplane.cpp
  Code/LB2.cpp
  Code/Lub.cpp
  Code/Model.cpp
  Code/Node.cpp
  Code/Point.cpp
  Code/project2.cpp
  Code/SLUB.cpp
  Code/solution.cpp
  Code/Tree.cpp
  Code/UB.cpp
)

add_executable(forget21 ${source_files})
target_link_libraries(forget21
  ${CPLEX_STUDIO}/cplex/lib/x86-64_linux/static_pic/libcplex.a
  ${CPLEX_STUDIO}/cplex/lib/x86-64_linux/static_pic/libilocplex.a
  ${CPLEX_STUDIO}/concert/lib/x86-64_linux/static_pic/libconcert.a
)
