cmake_minimum_required(VERSION 3.12)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR})
# Name your project and specify the required C++ standard
project(DiWaCAT_library)

# Specify the required C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Python and Pybind11
find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

# Define your extension
add_library(DiWaCAT_library MODULE FunctionsPythonBinding.cpp)

# Link your extension with Python and Pybind11
target_link_libraries(DiWaCAT_library PRIVATE pybind11::module Python3::Python)

# Optional: Specify the output directory for your .pyd file
set_target_properties(DiWaCAT_library PROPERTIES SUFFIX ".pyd")
set_target_properties(DiWaCAT_library PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
