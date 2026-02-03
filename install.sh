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
mkdir -p build && cd build &&
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake 2>/dev/null;  # Clean previous config
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=ON -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF .. && 
make -j8 && make install &&
cd .. &&

#compile the clust-mst
cd build &&
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake 2>/dev/null;  # Clean previous config
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF .. &&
make -j8 && make install &&
cd .. &&

#compile the clust-dbscan
cd build &&
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake 2>/dev/null;  # Clean previous config
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=ON .. &&
make -j8 && make install &&
cd ..

#compile the clust-leiden (DISABLED - removed leiden compilation)
# Leiden compilation has been removed from install.sh 
