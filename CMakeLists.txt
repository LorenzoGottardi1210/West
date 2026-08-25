#Same version of CMake required for Quantum ESPRESSO compilation 
cmake_minimum_required(VERSION 3.20, FATAL_ERROR)
project(WEST  VERSION 6.3.0 DESCRIPTION "WEST (Without Empty STates) is a massively parallel software for large scale electronic structure calculations within many-body perturbation theory "
        LANGUAGES C Fortran) 
#setting the enviromento to dump all the executables in the bin folder of the build directory of QE  

message(STATUS "\n=======================================================")
message(STATUS ">>> STARTING CONFIGURATION OF QUANTUM ESPRESSO <<<")
message(STATUS "=======================================================\n")
get_filename_component(QE_BIN_DIR "${CMAKE_SOURCE_DIR}/../bin" ABSOLUTE)
get_filename_component(UNIFIED_BIN_DIR "${CMAKE_SOURCE_DIR}/../bin" ABSOLUTE)

# setting path for the executables to be dumped in the bin folder of the build directory of QE
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UNIFIED_BIN_DIR} CACHE PATH "Unified bin directory" FORCE)

message(STATUS ">>> ALL EXECUTABLES WILL BE DUMPED TO: ${UNIFIED_BIN_DIR} <<<")


#sanity check before starting the compilation of WEST

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/qe_build)
set(Python3_FIND_VIRTUALENV FIRST) 

#Extracting compilation option from Quantum ESPRESSO and setting them for WEST compilation
get_directory_property(QE_FINAL_FORTRAN_FLAGS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. DEFINITION CMAKE_Fortran_FLAGS)
get_directory_property(QE_FINAL_C_FLAGS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. DEFINITION CMAKE_C_FLAGS)

if(QE_FINAL_FORTRAN_FLAGS)
    # setting compilation options for WEST to be the same as those used for Quantum ESPRESSO
    set(CMAKE_Fortran_FLAGS "${QE_FINAL_FORTRAN_FLAGS}")
    message(STATUS "WEST: Inherited FULL Fortran flags from QE: ${CMAKE_Fortran_FLAGS}")
endif()
if(QE_FINAL_C_FLAGS)
    set(CMAKE_C_FLAGS "${QE_FINAL_C_FLAGS}")
    message(STATUS "WEST: Inherited FULL C flags from QE: ${CMAKE_C_FLAGS}")
endif()

# using the target property to get the compilation options for WEST from Quantum ESPRESSO
get_target_property(QE_COMPILE_OPTS qe_modules COMPILE_OPTIONS)
if(NOT QE_COMPILE_OPTS STREQUAL "QE_COMPILE_OPTS-NOTFOUND")
    add_compile_options(${QE_COMPILE_OPTS})
endif()
#manual handling of CUDA and OpenACC flags for WEST compilation, as they are not automatically inherited from Quantum ESPRESSO
if(QE_ENABLE_CUDA AND CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    # searching for specific CUDA architecture flag set by the user for Quantum ESPRESSO compilation 
    if(DEFINED NVFORTRAN_CUDA_CC)
        set(GPU_ARCH_FLAG "-gpu=cc${NVFORTRAN_CUDA_CC}")
    else()
        # generic CUDA architecture flag for WEST compilation if not set by the user for Quantum ESPRESSO compilation 
        set(GPU_ARCH_FLAG "-Mcuda")
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda ${GPU_ARCH_FLAG}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -cuda ${GPU_ARCH_FLAG}")
    message(STATUS "WEST: Restored hidden NVHPC GPU flags (-cuda ${GPU_ARCH_FLAG})")
endif()

if(QE_ENABLE_OPENACC AND CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc")
    message(STATUS "WEST: Restored hidden NVHPC OpenACC flag (-acc)")
endif()

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
#The order here is important, as the libraries are needed for the modules to compile correctly
add_subdirectory(Libraries/Forpy)
add_subdirectory(Libraries/Base64)
add_subdirectory(Libraries/Json)
add_subdirectory(Modules)
add_subdirectory(Tools)
add_subdirectory(FFT_kernel) 
add_subdirectory(Hamiltonian_kernel)
add_subdirectory(Para_kernel)
add_subdirectory(Coulomb_kernel)
add_subdirectory(IO_kernel)
add_subdirectory(DFPT_kernel)

#lasts executables to be compiled, as they depend on all the libraries and modules above 
add_subdirectory(Wstat)
add_subdirectory(Wfreq)
add_subdirectory(Wbse)