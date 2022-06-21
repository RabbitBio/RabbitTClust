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

#compile the RabbitTClust
mkdir -p build && cd build &&
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=ON .. && 
make -j8 && make install &&
cd ../ &&

#mkdir -p build && cd build &&
cd build &&
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF .. &&
make -j8 && make install &&
cd ../ 
