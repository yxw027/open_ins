#cmake_minimum_required（VERSION 2.8）
project(INSCore)
include_directories(include)
include_directories(./)


set(CAMKE_BUILD_TYPE RELEASE)

file(GLOB core_src  *.c src/*.c )
file(GLOB core_hdr *.h include/*.h)
add_library(${PROJECT_NAME}  STATIC ${core_src} ${core_hdr})