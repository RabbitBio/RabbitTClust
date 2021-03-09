# RabbitClustering

## RabbitTClust
RabbitTClust is a clustering method based on Minimum Spanning Tree. The distance computing is based on Sketch.

Starting:
```bash
git clone --recursive git@github.com:RabbitBio/RabbitTClust.git
cd RabbitSketch
mkdir build
cd build 
cmake -DCXX_API=ON -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH

cd ../../
make 

./clust -l -t 4 -d 0.3 fileList >result

```

