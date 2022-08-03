cmake_minimum_required (VERSION 3.4)
project (kTests C CXX)
set(default_build_type "RELEASE")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC -lgomp -lrt -fopenmp -Ofast")

set(kProcessor "${PROJECT_SOURCE_DIR}/lib/kProcessor")

# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(${kProcessor})
set(kProcessor_INCLUDE_PATH "${kProcessor}/include")

# Add include directories.
include_directories("${kProcessor}/ThirdParty/MQF/include")
include_directories("${kProcessor}/ThirdParty/Blight")
include_directories (${kProcessor_INCLUDE_PATH})
include_directories("${kProcessor}/ThirdParty/caches/include")


# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------


add_executable (bins_indexing index_bins.cpp)
target_link_libraries (bins_indexing kProcessor z)
target_include_directories(bins_indexing INTERFACE ${kProcessor_INCLUDE_PATH})
