set -x

#make rabbitSketch library
cd RabbitSketch &&
mkdir -p build && cd build &&
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. .. &&
make && make install &&
cd ../../ &&

#make rabbitFX library
cd RabbitFX && 
mkdir -p build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make && make install && 
cd ../../ &&

#compile the RabbitTClust
mkdir -p build && cd build &&
cmake -DUSE_RABBITFX=ON .. && 
make && make install &&
cd ../ 
