
FILE(GLOB CPP_FILES *.cpp */*.cpp)

if(MSVC)
    set(CMAKE_CXX_FLAGS "-Wall -pedantic ${CMAKE_CXX_FLAGS}")
else(MSVC)
    set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic ${CMAKE_CXX_FLAGS}")
endif(MSVC)

add_library(Dimacs STATIC
    ${CPP_FILES}
    )

