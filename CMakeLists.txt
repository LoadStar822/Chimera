cmake_minimum_required(VERSION 3.10)
project(Chimera VERSION 1.6.3)
add_definitions(-DCHIMERA_VERSION="${PROJECT_VERSION}")
set(PROJECT_ROOT_DIR "${CMAKE_SOURCE_DIR}")
add_definitions(-DPROJECT_ROOT_DIR="${PROJECT_ROOT_DIR}")


set(CMAKE_BUILD_TYPE Release)

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

if( NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" )
    message( FATAL_ERROR
        "Compiler id '${CMAKE_CXX_COMPILER_ID}' is not supported, please \
        check the documentation." )
endif()



# 设置编译选项
include(CheckCXXCompilerFlag)

# 检测是否支持 AVX2 指令集
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

if(COMPILER_SUPPORTS_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2 -flto")
    message(STATUS "AVX2 instructions are supported and enabled.")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -flto")
    message(WARNING "AVX2 instructions are not supported. Proceeding without AVX2 optimization.")
endif()

# 设置SeqAn3库的路径
set(SEQAN3_DIR "${CMAKE_SOURCE_DIR}/libs/seqan3/share/cmake/seqan3")
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/libs/seqan3/share/cmake/seqan3")
set(SEQAN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/libs/seqan3/include")

# 查找 OpenSSL 包
find_package(OpenSSL REQUIRED)
find_package( BZip2 REQUIRED )
find_package( ZLIB REQUIRED )
find_package(Threads REQUIRED)
if (OpenSSL_FOUND)
    message(STATUS "Found OpenSSL: ${OPENSSL_VERSION}")
    message(STATUS "OpenSSL Include Dir: ${OPENSSL_INCLUDE_DIR}")
    message(STATUS "OpenSSL Libraries: ${OPENSSL_LIBRARIES}")
else()
    message(FATAL_ERROR "OpenSSL not found")
endif()
find_package (seqan3 REQUIRED)

# OpenMP（可选）
find_package(OpenMP)


# 包含目录
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs/
    ${PROJECT_SOURCE_DIR}/libs/klib
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
    ${PROJECT_SOURCE_DIR}/src/build/filter
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
list(FILTER SOURCES EXCLUDE REGEX "/tests/")
list(FILTER SOURCES EXCLUDE REGEX "/bench/")
list(FILTER SOURCES EXCLUDE REGEX "/packed/")
list(FILTER SOURCES EXCLUDE REGEX "/\\.worktrees/")
list(FILTER SOURCES EXCLUDE REGEX "/libs/minions-master/")
list(FILTER SOURCES EXCLUDE REGEX "/src/legacy/")
list(FILTER SOURCES EXCLUDE REGEX "/libs/strobemers/")
list(FILTER SOURCES EXCLUDE REGEX "/third_party/")

list(APPEND SOURCES "${PROJECT_SOURCE_DIR}/src/utils/FeatureHasher.cpp")
list(REMOVE_DUPLICATES SOURCES)

add_library(strobemers STATIC
    ${PROJECT_SOURCE_DIR}/libs/strobemers/randstrobes.cpp
)

target_include_directories(strobemers
    PUBLIC
    ${PROJECT_SOURCE_DIR}/libs/strobemers
)

# 递归地查找头文件
file(GLOB_RECURSE HEADERS
    "${PROJECT_SOURCE_DIR}/src/*.h"
    "${PROJECT_SOURCE_DIR}/src/*.hpp"
    "${PROJECT_SOURCE_DIR}/*.h"
    "${PROJECT_SOURCE_DIR}/*.hpp"
)

# 添加可执行文件
add_executable(Chimera ${SOURCES} "src/build/ChimeraBuild.cpp" "libs/xxhash.c" "libs/klib/kthread.c")


# 链接库和 OpenSSL
target_link_libraries(Chimera
    PRIVATE
    strobemers
    OpenSSL::SSL
    OpenSSL::Crypto
    pthread
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

target_compile_definitions(Chimera
    PRIVATE
    CHIMERA_HAS_STROBEMERS=1
)

# 设置包含目录
target_include_directories(Chimera
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/klib
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
    ${PROJECT_SOURCE_DIR}/src/build/filter
)




# IMCF 单元测试可执行文件
add_executable(imcf_filter_tests
    tests/imcf/interleaved_merged_cuckoo_filter_test.cpp
    src/utils/FeatureHasher.cpp
    libs/xxhash.c
    libs/klib/kthread.c
)

# 复用与主程序相同的依赖与包含路径
target_link_libraries(imcf_filter_tests
    PRIVATE
    strobemers
    OpenSSL::SSL
    OpenSSL::Crypto
    pthread
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

target_compile_definitions(imcf_filter_tests
    PRIVATE
    CHIMERA_HAS_STROBEMERS=1
)

if (OpenMP_CXX_FOUND)
    target_link_libraries(Chimera PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(imcf_filter_tests PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found: enabling parallel build.")
else()
    message(WARNING "OpenMP not found: building without parallelism.")
endif()

target_include_directories(imcf_filter_tests
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/build
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/klib
    ${PROJECT_SOURCE_DIR}/libs/minions-master/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
    ${PROJECT_SOURCE_DIR}/src/build/filter
)

add_executable(syncmer_tests
    tests/utils/syncmer_test.cpp
    src/utils/FeatureHasher.cpp
    libs/xxhash.c
)

target_link_libraries(syncmer_tests
    PRIVATE
    strobemers
    OpenSSL::SSL
    OpenSSL::Crypto
    pthread
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

target_compile_definitions(syncmer_tests
    PRIVATE
    CHIMERA_HAS_STROBEMERS=1
)

if (OpenMP_CXX_FOUND)
    target_link_libraries(syncmer_tests PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(syncmer_tests
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/minions-master/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
)

add_executable(unique_edge_guard_tests
    tests/utils/unique_edge_guard_test.cpp
)

target_link_libraries(unique_edge_guard_tests
    PRIVATE
    strobemers
    OpenSSL::SSL
    OpenSSL::Crypto
    pthread
    BZip2::BZip2
    ZLIB::ZLIB
    seqan3::seqan3
    Threads::Threads
)

target_include_directories(unique_edge_guard_tests
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/classify
    ${PROJECT_SOURCE_DIR}/libs
    ${PROJECT_SOURCE_DIR}/libs/minions-master/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/cereal/include
    ${PROJECT_SOURCE_DIR}/libs/seqan3/include/seqan3/submodules/sdsl-lite/include
)

# 忽略 register 关键字警告
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wno-register)
    add_compile_options(-Wno-sign-compare) # 忽略 signed/unsigned 比较警告
endif()


# 安装规则
install(TARGETS Chimera
    DESTINATION bin
)
