mkdir -p debug
cd debug
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug .. -DCMAKE_INSTALL_PREFIX=. -DNO_OPENMP=True && make install
