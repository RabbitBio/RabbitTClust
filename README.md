# RabbitTClust

RabbitTClust is a clustering method based on Minimum Spanning Tree. The distance computing is based on Sketch.

## Installation

### Dependancy
* c++14
* [zlib](https://zlib.net/)

Starting:
```bash
git clone --recursive git@github.com:RabbitBio/RabbitTClust.git
#make rabbitSketch library
cd RabbitSketch
mkdir build && cd build 
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
cd ../../

#make rabbitIO library
cd RabbitIO
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
cd ../../

mkdir build && cd build
cmake -DUSE_RABBITIO=ON ..
make && make install
cd ../

#The refList is the list path of the RefSeq genome files.
./clust -l -t 48 -i refList -o ref.out

#get the clustering result by inputing MST info.
./clust -f -d 0.01 -i refListMinHashGenomeInfo refListMinHashMSTInfo -o result.out

#get more help info.
./clust -h

```

## Experiment
The experiment idea and setting are in [EXPERIMENT](./experiment.md)
