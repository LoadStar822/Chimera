cmake_minimum_required(VERSION 3.10)

project(Chimera)

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

if( NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    message( FATAL_ERROR
        "Compiler id '${CMAKE_CXX_COMPILER_ID}' is not supported, please \
        check the documentation." )
endif()



# 设置编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-strict-aliasing -O3")

# 查找 OpenSSL 包
find_package(OpenSSL REQUIRED)

# 包含目录
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/lib/seqan3/include
)

# 递归地查找源文件，但排除 CMake 生成的文件
file(GLOB_RECURSE SOURCES
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.cc"
    "${PROJECT_SOURCE_DIR}/*.cpp"
    "${PROJECT_SOURCE_DIR}/*.cc"
)

# 移除 CMake 生成的文件
list(FILTER SOURCES EXCLUDE REGEX "CMakeFiles/")

# 递归地查找头文件
file(GLOB_RECURSE HEADERS
    "${PROJECT_SOURCE_DIR}/src/*.h"
    "${PROJECT_SOURCE_DIR}/src/*.hpp"
    "${PROJECT_SOURCE_DIR}/*.h"
    "${PROJECT_SOURCE_DIR}/*.hpp"
)

# 添加可执行文件
add_executable(Chimera ${SOURCES})

# 链接库和 OpenSSL
target_link_libraries(Chimera
    PRIVATE
    OpenSSL::SSL
    OpenSSL::Crypto
    pthread
)

# 设置包含目录
target_include_directories(Chimera
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/lib/seqan3/include
)

# 安装规则
install(TARGETS Chimera
    DESTINATION bin
)
