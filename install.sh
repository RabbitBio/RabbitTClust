set -euxo pipefail

# Always run from the directory containing this script (so paths are stable).
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Main RabbitTClust CMakeLists uses in-source style paths (RabbitSketch/build, ...).
# Build *out of tree* in these dirs so CMake never needs to write root-owned CMakeFiles/
# under the repo root (e.g. after a past `sudo cmake`).
: "${CMAKE_BUILD_TYPE:=Release}"
: "${MAKE_JOBS:=$(nproc 2>/dev/null || echo 8)}"

# Optional: remove stale in-source CMake junk if you own it (ignore errors).
rm -f CMakeCache.txt cmake_install.cmake Makefile 2>/dev/null || true
rm -rf CMakeFiles 2>/dev/null || true

# --- RabbitSketch ---
cd RabbitSketch
mkdir -p build && cd build
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. ..
make -j"${MAKE_JOBS}" && make install
cd "$SCRIPT_DIR"

# --- RabbitFX ---
cd RabbitFX
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j"${MAKE_JOBS}" && make install
cd "$SCRIPT_DIR"

# --- clust-greedy (out-of-tree; CMAKE_INSTALL_PREFIX .. = project root) ---
rm -rf cmake-build-greedy
mkdir -p cmake-build-greedy && cd cmake-build-greedy
cmake -DCMAKE_INSTALL_PREFIX=.. \
	-DUSE_RABBITFX=ON -DUSE_GREEDY=ON -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF ..
make -j"${MAKE_JOBS}" && make install
cd "$SCRIPT_DIR"

# --- clust-mst ---
rm -rf cmake-build-mst
mkdir -p cmake-build-mst && cd cmake-build-mst
cmake -DCMAKE_INSTALL_PREFIX=.. \
	-DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=OFF ..
make -j"${MAKE_JOBS}" && make install
cd "$SCRIPT_DIR"

# --- clust-dbscan ---
rm -rf cmake-build-dbscan
mkdir -p cmake-build-dbscan && cd cmake-build-dbscan
cmake -DCMAKE_INSTALL_PREFIX=.. \
	-DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=ON ..
make -j"${MAKE_JOBS}" && make install
cd "$SCRIPT_DIR"

echo "Done. Binaries installed under: $SCRIPT_DIR (see clust-greedy, clust-mst, clust-dbscan)"
