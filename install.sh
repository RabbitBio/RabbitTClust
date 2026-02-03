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

#compile the clust-dbscan
cd build &&
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=OFF -DUSE_DBSCAN=ON .. &&
make -j8 && make install &&
cd ../ &&

#compile the clust-leiden (requires igraph to be installed)
# Check if igraph is available (system-wide or in common locations)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "====================================================="
echo "Checking for igraph library..."
echo "====================================================="

# Priority 1: User-specified IGRAPH_ROOT environment variable
if [ -n "$IGRAPH_ROOT" ]; then
    echo "Checking user-specified IGRAPH_ROOT: $IGRAPH_ROOT"
    if [ -f "$IGRAPH_ROOT/lib64/pkgconfig/igraph.pc" ]; then
        export PKG_CONFIG_PATH="$IGRAPH_ROOT/lib64/pkgconfig:$PKG_CONFIG_PATH"
        export LD_LIBRARY_PATH="$IGRAPH_ROOT/lib64:$LD_LIBRARY_PATH"
        echo "✓ Found igraph in: $IGRAPH_ROOT/lib64"
    elif [ -f "$IGRAPH_ROOT/lib/pkgconfig/igraph.pc" ]; then
        export PKG_CONFIG_PATH="$IGRAPH_ROOT/lib/pkgconfig:$PKG_CONFIG_PATH"
        export LD_LIBRARY_PATH="$IGRAPH_ROOT/lib:$LD_LIBRARY_PATH"
        echo "✓ Found igraph in: $IGRAPH_ROOT/lib"
    else
        echo "✗ igraph.pc not found in IGRAPH_ROOT: $IGRAPH_ROOT"
    fi
fi

# Priority 2: Check if already in PKG_CONFIG_PATH (e.g., system-wide installation)
if pkg-config --exists igraph 2>/dev/null; then
    IGRAPH_LOCATION=$(pkg-config --variable=prefix igraph 2>/dev/null)
    echo "✓ igraph found via pkg-config: $IGRAPH_LOCATION"
else
    # Priority 3: Search common local installation paths
    echo "Searching common installation paths..."
    SEARCH_PATHS=(
        "$IGRAPH_ROOT"
        "$HOME/leiden/igraph-install"
        "$SCRIPT_DIR/../../../igraph-install"
        "$SCRIPT_DIR/../../igraph-install"
        "$SCRIPT_DIR/../igraph-install"
        "$HOME/igraph-install"
        "$HOME/local"
        "$HOME/.local"
        "/usr/local"
        "/opt/igraph"
        "/opt/local"
    )
    
    for IGRAPH_DIR in "${SEARCH_PATHS[@]}"; do
        [ -z "$IGRAPH_DIR" ] && continue
        if [ -f "$IGRAPH_DIR/lib64/pkgconfig/igraph.pc" ]; then
            export PKG_CONFIG_PATH="$IGRAPH_DIR/lib64/pkgconfig:$PKG_CONFIG_PATH"
            export LD_LIBRARY_PATH="$IGRAPH_DIR/lib64:$LD_LIBRARY_PATH"
            echo "✓ Found igraph in: $IGRAPH_DIR/lib64"
            break
        elif [ -f "$IGRAPH_DIR/lib/pkgconfig/igraph.pc" ]; then
            export PKG_CONFIG_PATH="$IGRAPH_DIR/lib/pkgconfig:$PKG_CONFIG_PATH"
            export LD_LIBRARY_PATH="$IGRAPH_DIR/lib:$LD_LIBRARY_PATH"
            echo "✓ Found igraph in: $IGRAPH_DIR/lib"
            break
        fi
    done
fi

# Final check: verify igraph is accessible
if pkg-config --exists igraph; then
    IGRAPH_VERSION=$(pkg-config --modversion igraph 2>/dev/null || echo "unknown")
    IGRAPH_PREFIX=$(pkg-config --variable=prefix igraph 2>/dev/null || echo "unknown")
    echo "====================================================="
    echo "✓ igraph detected"
    echo "  Version: $IGRAPH_VERSION"
    echo "  Location: $IGRAPH_PREFIX"
    echo "Compiling clust-leiden..."
    echo "====================================================="
    cd build &&
    cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF -DUSE_LEIDEN=ON .. &&
    make -j8 && make install &&
    cd ../
    echo ""
    echo "====================================================="
    echo "✓ clust-leiden compiled successfully"
    echo "====================================================="
else
    echo "====================================================="
    echo "✗ WARNING: igraph library not detected"
    echo ""
    echo "clust-leiden will NOT be compiled."
    echo "clust-greedy and clust-mst are still available."
    echo ""
    echo "To enable clust-leiden:"
    echo ""
    echo "1. If igraph is already installed elsewhere:"
    echo "   Edit install.sh, line ~62"
    echo "   Change \"\\\$HOME/leiden/igraph-install\""
    echo "   to your igraph installation path"
    echo ""
    echo "2. If igraph is not installed:"
    echo "   sudo apt-get install libigraph-dev  # Ubuntu/Debian"
    echo "   sudo yum install igraph-devel       # CentOS/RHEL"
    echo ""
    echo "====================================================="
fi 
