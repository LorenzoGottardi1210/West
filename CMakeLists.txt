#Same version of CMake required for Quantum ESPRESSO compilation 
cmake_minimum_required(VERSION 3.20, FATAL_ERROR)
project(WEST 
        VERSION 6.3.0 
        DESCRIPTION "WEST (Without Empty STates) is a massively parallel software
                            for large scale electronic structure calculations within many-body perturbation theory "
        LANGUAGES C Fortran) 
#setting the enviromento to dump all the executables in the bin folder of the build directory of QE  

message(STATUS "\n=======================================================")
message(STATUS ">>> STARTING CONFIGURATION OF QUANTUM ESPRESSO <<<")
message(STATUS "=======================================================\n")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#sanity check before starting the compilation of WEST

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/qe_build)
set(Python3_FIND_VIRTUALENV FIRST) 

find_package(Python3  COMPONENTS Interpreter Development)


message(STATUS "=======================================================")
message(STATUS "Python found for WEST:")
message(STATUS "  - Executable: ${Python3_EXECUTABLE}")
message(STATUS "  - Library (Din. link): ${Python3_LIBRARIES}")
message(STATUS "  - Header (Include): ${Python3_INCLUDE_DIRS}")
message(STATUS "=======================================================") 
message(STATUS "Configuration of QE terminated, proceeding to WEST <<<")



if(NOT Python3_Development_FOUND)
    message(WARNING "WEST requires Python3 development libraries to be installed. Please install the appropriate package for your system (e.g., python3-dev or python3-devel) and re-run CMake.")
endif() 
#actually compiling WEST 
message(STATUS ">>> COMPILING WEST  <<<")
#First the subdirectory Pytools is compiled, then the libraries Forpy, Base64 and Json are compiled. They doesn't have any dependencies between them, so the order of compilation doesn't matter.
add_subdirectory(Pytools)
add_subdirectory(Libraries/Forpy)
add_subdirectory(Libraries/Base64)
add_subdirectory(Libraries/Json)