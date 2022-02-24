set -x
cd RabbitSketch &&
mkdir build && cd build &&
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. .. &&
make && make install &&
#export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH &&
cd ../../ &&

#make rabbitFX library
cd RabbitFX && 
mkdir build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make && make install && 
#export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH &&
cd ../../ &&

mkdir build && cd build &&
cmake -DUSE_RABBITFX=ON .. && 
make && make install &&
cd ../ 
