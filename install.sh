set -x

#make rabbitSketch library
cd RabbitSketch &&
mkdir -p build && cd build &&
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install &&
cd ../../ &&

#make rabbitFX library
cd RabbitFX && 
mkdir -p build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install && 
cd ../../ &&

#compile the clust-greedy
mkdir -p build &&
rm -rf build/CMakeCache.txt build/CMakeFiles build/Makefile build/cmake_install.cmake 2>/dev/null;  # Clean build dir config
cmake -S . -B build -DUSE_RABBITFX=ON -DUSE_GREEDY=ON -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF &&
cmake --build build -j8 && cmake --install build &&

#compile the clust-mst
rm -rf build/CMakeCache.txt build/CMakeFiles build/Makefile build/cmake_install.cmake 2>/dev/null;  # Clean build dir config
cmake -S . -B build -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF &&
cmake --build build -j8 && cmake --install build &&

#compile the clust-dbscan
rm -rf build/CMakeCache.txt build/CMakeFiles build/Makefile build/cmake_install.cmake 2>/dev/null;  # Clean build dir config
cmake -S . -B build -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=ON &&
cmake --build build -j8 && cmake --install build

#compile the clust-leiden (DISABLED - removed leiden compilation)
# Leiden compilation has been removed from install.sh 
