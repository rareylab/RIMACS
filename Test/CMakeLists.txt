function(build_test TEST_NAME TEST_FILE FLAGS)
  set(TEST_EXECUTABLE "${TEST_NAME}_tester")
  add_executable(${TEST_EXECUTABLE} main.cpp ${TEST_FILE})

  set_property(TARGET ${TEST_EXECUTABLE} PROPERTY CXX_STANDARD 14)
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES GNU AND ${RIMACS_WITH_COV})
    set_target_properties("${TEST_EXECUTABLE}" PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -coverage"
                          LINK_FLAGS "-coverage")
  endif()

  if(MSVC)
    set_target_properties("${TEST_EXECUTABLE}" PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
    target_link_libraries (${TEST_EXECUTABLE} RIMACSStatic Dimacs)
  else(MSVC)
    target_link_libraries (${TEST_EXECUTABLE} RIMACS Dimacs)
  endif(MSVC)
    target_link_libraries (${TEST_EXECUTABLE} gtest_main)

  add_test(NAME ${TEST_NAME} COMMAND ${TEST_EXECUTABLE} ${FLAGS})
  add_dependencies(build_and_test ${TEST_NAME}_tester)
endfunction()

# add subdirectory for simple graph representation
add_subdirectory(Dimacs)

build_test(SubgraphMapping SubgraphMappingTester.cpp "")
build_test(VertexMapping VertexMappingTester.cpp "")
build_test(MCS MCSTester.cpp "")
build_test(MCES MCESTester.cpp "")
build_test(CachingFunctor FunctorTester.cpp "")
build_test(MCSResult MCSResultTester.cpp "")
build_test(Example RIMACSexampleTester.cpp "")
build_test(EquivalenceClass EquivalenceClassTester.cpp "")
build_test(NodeOrder NodeOrderTester.cpp "")
build_test(MinimumConstraints MinimumResultConstraintTester.cpp "")

# copy all test dimacs files to execution folder once at cmake configuration time
if(MSVC)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Examples/
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/Examples)
else(MSVC)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Examples/
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Examples)
endif(MSVC)
