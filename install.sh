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
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=ON -DUSE_LEIDEN=OFF .. && 
make -j8 && make install &&
cd ../ &&

#compile the clust-mst
cd build &&
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF .. &&
make -j8 && make install &&
cd ../ &&

#compile the clust-leiden (requires igraph to be installed)
# Check if igraph is available (system-wide or in common locations)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Try to find igraph in common local installation paths
# Check for igraph.pc file to ensure igraph is actually there
for IGRAPH_DIR in \
    "$SCRIPT_DIR/../../../igraph-install" \
    "$SCRIPT_DIR/../../igraph-install" \
    "$HOME/local" \
    "$HOME/.local" \
    "/usr/local" \
    "/opt/igraph"
do
    if [ -f "$IGRAPH_DIR/lib64/pkgconfig/igraph.pc" ]; then
        export PKG_CONFIG_PATH="$IGRAPH_DIR/lib64/pkgconfig:$PKG_CONFIG_PATH"
        export LD_LIBRARY_PATH="$IGRAPH_DIR/lib64:$LD_LIBRARY_PATH"
        echo "Found igraph in: $IGRAPH_DIR"
        break
    elif [ -f "$IGRAPH_DIR/lib/pkgconfig/igraph.pc" ]; then
        export PKG_CONFIG_PATH="$IGRAPH_DIR/lib/pkgconfig:$PKG_CONFIG_PATH"
        export LD_LIBRARY_PATH="$IGRAPH_DIR/lib:$LD_LIBRARY_PATH"
        echo "Found igraph in: $IGRAPH_DIR"
        break
    fi
done

# Check if igraph is available via pkg-config
if pkg-config --exists igraph; then
    IGRAPH_VERSION=$(pkg-config --modversion igraph 2>/dev/null || echo "unknown")
    echo "====================================================="
    echo "igraph detected (version: $IGRAPH_VERSION)"
    echo "Compiling clust-leiden..."
    echo "====================================================="
    cd build &&
    cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=ON .. &&
    make -j8 && make install &&
    cd ../
    echo "✓ clust-leiden compiled successfully"
else
    echo "====================================================="
    echo "WARNING: igraph library not detected"
    echo "Skipping clust-leiden compilation"
    echo "Note: clust-greedy and clust-mst are available"
    echo "====================================================="
fi 
