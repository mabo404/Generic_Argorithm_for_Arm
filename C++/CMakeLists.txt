cmake_minimum_required(VERSION 3.15)
project(nsga_arm CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## --- pagmo をヘッダ/ライブラリ直接探索に切り替え ---
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
endif()
find_path(PAGMO_INCLUDE_DIR
  NAMES pagmo/algorithms/nsga2.hpp
  HINTS "$ENV{CONDA_PREFIX}/include" "${CMAKE_PREFIX_PATH}/include")
find_library(PAGMO_LIBRARY
  NAMES pagmo
  HINTS "$ENV{CONDA_PREFIX}/lib" "${CMAKE_PREFIX_PATH}/lib")
if(NOT PAGMO_INCLUDE_DIR OR NOT PAGMO_LIBRARY)
  message(FATAL_ERROR "pagmo headers or library not found. Set CONDA_PREFIX or CMAKE_PREFIX_PATH.")
endif()

add_executable(nsga main.cpp)
target_include_directories(nsga PRIVATE ${PAGMO_INCLUDE_DIR})
target_link_libraries(nsga PRIVATE ${PAGMO_LIBRARY})
# 実行時に conda の .so を見つけやすくする（任意）
set_target_properties(nsga PROPERTIES BUILD_RPATH "$ENV{CONDA_PREFIX}/lib")
