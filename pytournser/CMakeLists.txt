cmake_minimum_required(VERSION 2.8.12)
project(pytournser)
add_subdirectory(pybind11)
pybind11_add_module(pytournser tournser_bindings.cpp)
