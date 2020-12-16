
FILE(GLOB CPP_FILES *.cpp */*.cpp)

if(MSVC)
    set(CMAKE_CXX_FLAGS "-Wall -pedantic ${CMAKE_CXX_FLAGS}")
else(MSVC)
    set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic ${CMAKE_CXX_FLAGS}")
endif(MSVC)

function(build_library LIBRARY_NAME SPECIFIER)
  add_library(${LIBRARY_NAME} ${SPECIFIER}
      ${CPP_FILES}
      )
  set_target_properties("${LIBRARY_NAME}" PROPERTIES POSITION_INDEPENDENT_CODE True)
  install(TARGETS ${LIBRARY_NAME}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

endfunction()

build_library("RIMACS" "SHARED")
build_library("RIMACSStatic" "STATIC")
# Handle install of header files
install(DIRECTORY "${INCLUDE_DIR}/RIMACS" DESTINATION include FILES_MATCHING PATTERN "*.hpp")

if((${CMAKE_CXX_COMPILER_ID} MATCHES GNU) AND ${RIMACS_WITH_COV})
  set_target_properties("RIMACS" PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -coverage"
                        LINK_FLAGS "-coverage")
endif()

