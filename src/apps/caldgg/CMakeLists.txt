
cmake_minimum_required (VERSION 3.9)

add_executable(caldgg main.cpp calDgg.cpp fileio.cpp testfun.h )
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(caldgg PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(caldgg PRIVATE dglib)
target_compile_features(caldgg PRIVATE cxx_std_11)





