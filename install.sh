set -x
cd RabbitSketch &&
mkdir build && cd build &&
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. .. &&
make && make install &&
cd ../../ &&

#make rabbitFX library
cd RabbitFX && 
mkdir build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make && make install && 
cd ../../ &&

mkdir build && cd build &&
cmake -DUSE_RABBITFX=ON .. && 
make && make install &&
cd ../ 
