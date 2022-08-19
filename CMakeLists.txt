cmake_minimum_required(VERSION 2.8.12)
project(pytourn)
add_subdirectory(./pytournser/pybind11)
pybind11_add_module(pytourn ./pytournser/tournser_bindings.cpp)
