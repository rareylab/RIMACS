cmake_minimum_required (VERSION 3.1)
project (RIMACS VERSION 1.0)

# enable the CMake testing
enable_testing()
option(run_tests "Run unit tests automatically" ON)

# perform coverage report
option(RIMACS_WITH_COV "Enable coverage report analysis" OFF)

# set the include directories for boost and Eigen
set(BOOST_INCLUDE_DIR "" CACHE STRING "Include dir where to find boost")
set(EIGEN_INCLUDE_DIR "" CACHE STRING "Include dir where to find Eigen")

set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")

# currently needed by google test on macOS
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -ggdb")

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()
# Add make target for test compilation and execution
# add -V if you like to see the full GTEST output
add_custom_target(build_and_test ${CMAKE_CTEST_COMMAND})

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

if (NOT (EIGEN_INCLUDE_DIR STREQUAL ""))
  include_directories(BEFORE SYSTEM ${EIGEN_INCLUDE_DIR})
  message("Included Eigen ${EIGEN_INCLUDE_DIR}")
else()
  find_package (Eigen3 3.3 REQUIRED NO_MODULE)
endif()

if (NOT (BOOST_INCLUDE_DIR STREQUAL ""))
  include_directories(BEFORE SYSTEM ${BOOST_INCLUDE_DIR})
  message("Included BOOST ${BOOST_INCLUDE_DIR}")
else()
  find_package(Boost 1.60.0 REQUIRED)
endif()

include_directories(${INCLUDE_DIR})
add_subdirectory("${SOURCE_DIR}/RIMACS")
add_subdirectory("${CMAKE_SOURCE_DIR}/Test")

# generate meaningful coverage report. Requires the lcov tool
if((${CMAKE_CXX_COMPILER_ID} MATCHES GNU)  AND ${RIMACS_WITH_COV})
  message(STATUS "Enabled coverage builds")

  add_custom_target(extract_coverage COMMAND /usr/bin/lcov --quiet --capture --output-file allCoverage.info --directory src/RIMACS --directory Test/
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" COMMENT "Extract all coverage information" DEPENDS "build_and_test"
                    USES_TERMINAL VERBATIM BYPRODUCTS allCoverage.info)
  add_custom_target(filter_coverage COMMAND lcov --quiet --extract allCoverage.info '*src/RIMACS*' '*include/RIMACS*' --output-file RIMACScoverage.info
                    DEPENDS "extract_coverage" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
                    COMMENT "Filter RIMACS coverage" BYPRODUCTS RIMACScoverage.info)
  add_custom_target(generate_report COMMAND genhtml --demangle-cpp RIMACScoverage.info --output-directory ${CMAKE_SOURCE_DIR}/RIMACS-coverage-report
                    DEPENDS "filter_coverage" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" COMMENT "Generate HTML report in SRC directory"
                    BYPRODUCTS ${CMAKE_SOURCE_DIR}/RIMACS-coverage-report)
  add_custom_target(clear_counters COMMAND lcov --zerocounters --directory . DEPENDS "extract_coverage"
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" COMMENT "Cleanup call counters")
  add_custom_target(coverage_report COMMAND echo Please run: \"sensible-browser ${CMAKE_SOURCE_DIR}/RIMACS-coverage-report/index.html\" to inspect the coverage report
                    DEPENDS "generate_report" "clear_counters" COMMENT "Execute whole report chain")
endif()
