# RabbitClustering

## RabbitTClust
RabbitTClust is a clustering method based on Minimum Spanning Tree. The distance computing is based on Sketch.

Starting:
```bash
git clone --recursive git@github.com:RabbitBio/RabbitTClust.git
cd RabbitSketch
mkdir build
cd build 
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH

cd ../../
make 
#The fileList is the list path of the genome files.
./clust -l -t 4 -d 0.3 fileList >result

#get the clustering result by inputing MST info.
./clust -f -d 0.01 fileListGenomeInfo fileListMSTInfo

#get more help info.
./clust -h

```

## Experiment
The experiment idea and setting are in [EXPERIMENT](./experiment.md)
