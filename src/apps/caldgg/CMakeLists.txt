
cmake_minimum_required (VERSION 3.9)

add_executable(caldgg main.cpp calDgg.cpp fileio.cpp testfun.h )

target_link_libraries(caldgg PRIVATE dglib)
target_compile_features(caldgg PRIVATE cxx_std_11)





