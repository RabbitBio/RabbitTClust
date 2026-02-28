cmake_minimum_required(VERSION 3.0)

project(rabbitTClust)

set(CMAKE_INSTALL_PREFIX ..)
option(USE_RABBITFX "parse input fasta file with RabbitFX" OFF)
option(USE_DEBUG "print the debug information" ON)
option(USE_Timer "print the time information" ON)
option(USE_GREEDY "use greedy incremental cluster" ON)
# option(USE_LEIDEN "use Leiden community detection cluster" OFF)  # DISABLED
option(USE_DBSCAN "use DBSCAN density-based clustering" OFF)
option(USE_MPI "build MPI version of clust-mst (MinHash and KSSD)" ON)


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

# if(USE_LEIDEN)  # DISABLED
# 	set(CMAKE_CXX_FLAGS "-D LEIDEN_CLUST ${CMAKE_CXX_FLAGS}")
# 	
# 	# Find igraph library (only dependency needed)
# 	find_package(PkgConfig REQUIRED)
# 	pkg_check_modules(IGRAPH REQUIRED igraph)
# 	
# 	# Add igraph include directories
# 	include_directories(${IGRAPH_INCLUDE_DIRS})
# 	link_directories(${IGRAPH_LIBRARY_DIRS})
# endif()

if(USE_DBSCAN)
	set(CMAKE_CXX_FLAGS "-D DBSCAN_CLUST ${CMAKE_CXX_FLAGS}")
endif()

if(USE_MPI)
	find_package(MPI REQUIRED)
	add_definitions(-DUSE_MPI)
	include_directories(${MPI_CXX_INCLUDE_PATH})
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

# Add igraph include directories for Leiden (DISABLED)
# if(USE_LEIDEN)
# 	include_directories(${IGRAPH_INCLUDE_DIRS})
# 	link_directories(${IGRAPH_LIBRARY_DIRS})
# endif()

##include_directories(src RabbitSketch/build/include RabbitFX/build/include)
#include_directories(src RabbitSketch/build/include)
#
##link_directories(RabbitSketch/build/lib RabbitFX/build/io)
#link_directories(RabbitSketch/build/lib)

aux_source_directory(src DIR_SRCS)

# Exclude leiden.cpp from compilation (always exclude since Leiden is disabled)
list(REMOVE_ITEM DIR_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/leiden.cpp")
# Also try relative path
list(REMOVE_ITEM DIR_SRCS "src/leiden.cpp")

if(USE_GREEDY)
add_executable(clust-greedy ${DIR_SRCS})
if(USE_RABBITFX)
	target_link_libraries(clust-greedy rabbitsketch_static rabbitfx z)
else()
	target_link_libraries(clust-greedy rabbitsketch_static z)
endif()
if(USE_MPI)
	target_link_libraries(clust-greedy ${MPI_CXX_LIBRARIES})
endif()

# elseif(USE_LEIDEN)  # DISABLED
# add_executable(clust-leiden ${DIR_SRCS})
# if(USE_RABBITFX)
# 	target_link_libraries(clust-leiden rabbitsketch_static rabbitfx z ${IGRAPH_LIBRARIES} openblas gfortran)
# else()
# 	target_link_libraries(clust-leiden rabbitsketch_static z ${IGRAPH_LIBRARIES} openblas gfortran)
# endif()

elseif(USE_DBSCAN)
add_executable(clust-dbscan ${DIR_SRCS})
if(USE_RABBITFX)
	target_link_libraries(clust-dbscan rabbitsketch_static rabbitfx z)
else()
	target_link_libraries(clust-dbscan rabbitsketch_static z)
endif()
if(USE_MPI)
	target_link_libraries(clust-dbscan ${MPI_CXX_LIBRARIES})
endif()

else()
add_executable(clust-mst ${DIR_SRCS})
if(USE_RABBITFX)
	target_link_libraries(clust-mst rabbitsketch_static rabbitfx z)
else()
	target_link_libraries(clust-mst rabbitsketch_static z)
endif()
if(USE_MPI)
	target_link_libraries(clust-mst ${MPI_CXX_LIBRARIES})
endif()
endif()



if(USE_GREEDY)
install(TARGETS clust-greedy DESTINATION ./)
# elseif(USE_LEIDEN)  # DISABLED
# install(TARGETS clust-leiden DESTINATION ./)
elseif(USE_DBSCAN)
install(TARGETS clust-dbscan DESTINATION ./)
else()
install(TARGETS clust-mst DESTINATION ./)
endif()

