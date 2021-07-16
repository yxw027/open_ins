cmake_minimum_required(VERSION 3.16)
project(INS)

add_definitions(-DBUG)

add_subdirectory(Core)

include_directories(cJSON)
include_directories(Core)
include_directories(./)


file(GLOB HEADERFS  configManger.h util.h)
file(GLOB SOURCEFS  configManger.c util.c)
add_executable(${PROJECT_NAME} main.c ${HEADERFS} ${SOURCEFS})
#target_include_directories(${PROJECT_NAME} PRIVATE  cJSON)
link_directories(cJSON)
target_link_libraries(${PROJECT_NAME}  PRIVATE  cjson.lib INSCore)