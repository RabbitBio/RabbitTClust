cmake_minimum_required(VERSION 3.0)

project(rabbitTClust)

set(CMAKE_INSTALL_PREFIX ..)
option(USE_RABBITFX "parse input fasta file with RabbitFX" OFF)
option(USE_DEBUG "print the debug information" ON)
option(USE_Timer "print the time information" ON)
option(USE_GREEDY "use greedy incremental cluster" ON)



find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
#message("openmp found")

set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(EXECUTABLE_OUTPUT_PATH .)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "-g -O3 -D THREADPOOL_MINHASH -D DEBUG -D Timer ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "-g -O3 -Wno-format -march=native -Wno-unused-result${CMAKE_CXX_FLAGS}")

if(USE_DEBUG)
	set(CMAKE_CXX_FLAGS "-D DEBUG ${CMAKE_CXX_FLAGS}")
endif()

if(USE_Timer)
	set(CMAKE_CXX_FLAGS "-D Timer ${CMAKE_CXX_FLAGS}")
endif()

if(USE_GREEDY)
	set(CMAKE_CXX_FLAGS "-D GREEDY_CLUST ${CMAKE_CXX_FLAGS}")
endif()

if(USE_RABBITFX)
	set(CMAKE_CXX_FLAGS "-D RABBIT_FX ${CMAKE_CXX_FLAGS}")
	include_directories(src RabbitSketch/build/include RabbitFX/build/include)
	link_directories(RabbitSketch/build/lib RabbitFX/build/lib)
else()
	set(CMAKE_CXX_FLAGS "-D THREADPOOL_MINHASH ${CMAKE_CXX_FLAGS}")
	include_directories(src RabbitSketch/build/include)
	link_directories(RabbitSketch/build/lib)
endif()
	
##include_directories(src RabbitSketch/build/include RabbitFX/build/include)
#include_directories(src RabbitSketch/build/include)
#
##link_directories(RabbitSketch/build/lib RabbitFX/build/io)
#link_directories(RabbitSketch/build/lib)

aux_source_directory(src DIR_SRCS)

if(USE_GREEDY)
add_executable(clust-greedy ${DIR_SRCS})
if(USE_RABBITFX)
	target_link_libraries(clust-greedy rabbitsketch_static rabbitfx z)
else()
	target_link_libraries(clust-greedy rabbitsketch_static z)
endif()

else()
add_executable(clust-mst ${DIR_SRCS})
if(USE_RABBITFX)
	target_link_libraries(clust-mst rabbitsketch_static rabbitfx z)
else()
	target_link_libraries(clust-mst rabbitsketch_static z)
endif()
endif()



if(USE_GREEDY)
install(TARGETS clust-greedy DESTINATION ./)
else()
install(TARGETS clust-mst DESTINATION ./)
endif()

