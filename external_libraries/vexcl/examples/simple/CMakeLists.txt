cmake_minimum_required(VERSION 3.1)
project(hello)

find_package(VexCL)

vexcl_add_executables(hello hello.cpp)
