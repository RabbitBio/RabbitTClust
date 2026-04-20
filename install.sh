set -x

# --- Detect MPI (optional) ---
# USE_MPI defaults to ON in CMakeLists.txt, but the installer must work on
# machines without MPI too. If mpicxx is reachable, enable MPI and point
# CMake at the discovered wrappers; otherwise disable it so the build still
# produces working clust-greedy / clust-mst / clust-dbscan binaries.
MPI_CXX_BIN=$(command -v mpicxx 2>/dev/null)
MPI_C_BIN=$(command -v mpicc 2>/dev/null)
if [ -n "$MPI_CXX_BIN" ] && [ -n "$MPI_C_BIN" ]; then
    MPI_FLAGS=(-DUSE_MPI=ON "-DMPI_CXX_COMPILER=$MPI_CXX_BIN" "-DMPI_C_COMPILER=$MPI_C_BIN")
    echo "[install.sh] MPI detected: $MPI_CXX_BIN / $MPI_C_BIN"
else
    MPI_FLAGS=(-DUSE_MPI=OFF)
    echo "[install.sh] MPI not found on PATH -- building without MPI support"
fi

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

#compile the clust-greedy (in-source build)
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake 2>/dev/null  # Clean previous config
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=ON -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF "${MPI_FLAGS[@]}" . &&
make -j8 && make install &&

#compile the clust-mst
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake 2>/dev/null  # Clean previous config
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF "${MPI_FLAGS[@]}" . &&
make -j8 && make install &&

#compile the clust-dbscan
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake 2>/dev/null  # Clean previous config
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=ON "${MPI_FLAGS[@]}" . &&
make -j8 && make install

#compile the clust-leiden (DISABLED - removed leiden compilation)
# Leiden compilation has been removed from install.sh 
