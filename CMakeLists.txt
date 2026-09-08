#Same version of CMake required for Quantum ESPRESSO compilation 
cmake_minimum_required(VERSION 3.20, FATAL_ERROR)
project(WEST  VERSION 6.4.0 DESCRIPTION "WEST (Without Empty STates) is a massively parallel software for large scale electronic structure calculations within many-body perturbation theory "
        LANGUAGES C Fortran) 

# ==============================================================================
# STEP 1: enviroment to dump all the executables in the bin folder of the build directory of Q
# ==============================================================================
message(STATUS "\n=======================================================")
message(STATUS ">>> STEP 1: configuration of the bin directory <<<")
message(STATUS "=======================================================\n")

get_filename_component(QE_BIN_DIR "${CMAKE_SOURCE_DIR}/../bin" ABSOLUTE)
get_filename_component(UNIFIED_BIN_DIR "${CMAKE_SOURCE_DIR}/../bin" ABSOLUTE)

# setting path for the executables to be dumped in the bin folder of the build directory of QE
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UNIFIED_BIN_DIR} CACHE PATH "Unified bin directory" FORCE)
message(STATUS ">>> ALL EXECUTABLES WILL BE DUMPED TO: ${UNIFIED_BIN_DIR} <<<")

# ==============================================================================
# STEP 2: searching for QE flags to match them perfectly 1:1 with WEST
# ==============================================================================
message(STATUS "\n=======================================================")
message(STATUS ">>> STEP 2: searching for QE flags  <<<")
message(STATUS "=======================================================\n")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/qe_build)

#sniffing the flags of QE to match them perfectly 1:1 with WEST 
get_directory_property(QE_BASE_F_FLAGS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. DEFINITION CMAKE_Fortran_FLAGS)
get_directory_property(QE_BASE_C_FLAGS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. DEFINITION CMAKE_C_FLAGS)
string(TOUPPER "${CMAKE_BUILD_TYPE}" BTYPE)
get_directory_property(QE_OPT_F_FLAGS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. DEFINITION CMAKE_Fortran_FLAGS_${BTYPE})
get_directory_property(QE_OPT_C_FLAGS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. DEFINITION CMAKE_C_FLAGS_${BTYPE})

# setting the flags for WEST to match perfectly 1:1 with QE 
set(CMAKE_Fortran_FLAGS "${QE_BASE_F_FLAGS} ${QE_OPT_F_FLAGS}")
set(CMAKE_C_FLAGS "${QE_BASE_C_FLAGS} ${QE_OPT_C_FLAGS}")

# printing the flags sniffed from QE for sanity check 
message(STATUS "WEST: Flag C found from  QE       -> ${CMAKE_C_FLAGS}")
message(STATUS "WEST: Flag Fortran found from QE  -> ${CMAKE_Fortran_FLAGS}")

# 2. MODERN CMAKE SNIFFING (capture target properties from QE)
get_target_property(QE_COMPILE_OPTS qe_modules COMPILE_OPTIONS)
if(QE_COMPILE_OPTS AND NOT QE_COMPILE_OPTS STREQUAL "QE_COMPILE_OPTS-NOTFOUND")
    add_compile_options(${QE_COMPILE_OPTS})
    message(STATUS "WEST: Target Properties found from QE   -> ${QE_COMPILE_OPTS}")
endif()

# 3. dynamical translation of the GPU architectures for NVHPC/PGI compilers
if(QE_GPU AND CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(WEST_GPU_FLAGS "")
    
    if(QE_GPU MATCHES "openacc")
        set(WEST_GPU_FLAGS "${WEST_GPU_FLAGS} -acc")
    endif()
    
    if(QE_GPU MATCHES "cuda")
        set(WEST_GPU_FLAGS "${WEST_GPU_FLAGS} -cuda")
        
        if(QE_GPU_ARCHS)
            string(REPLACE ";" "," NV_ARCH_LIST "${QE_GPU_ARCHS}")
            string(REPLACE "sm_" "cc" NV_ARCH_FINAL "${NV_ARCH_LIST}")
            set(WEST_GPU_FLAGS "${WEST_GPU_FLAGS} -gpu=${NV_ARCH_FINAL}")
        else()
            set(WEST_GPU_FLAGS "${WEST_GPU_FLAGS} -Mcuda")
        endif()
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${WEST_GPU_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${WEST_GPU_FLAGS}")
    add_link_options(${WEST_GPU_FLAGS})
    
    message(STATUS "WEST: Traduzione hardware NVHPC   -> ${WEST_GPU_FLAGS}")
endif()

# ==============================================================================
# STEP 3: Python interpreter and development libraries for WEST 
# ==============================================================================
message(STATUS "\n=======================================================")
message(STATUS ">>> STEP 3: configuration of Python for WEST <<<")
message(STATUS "=======================================================\n")
set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3  COMPONENTS Interpreter Development)
message(STATUS "Python Executable : ${Python3_EXECUTABLE}")
message(STATUS "Python Library    : ${Python3_LIBRARIES}")
if(NOT Python3_Development_FOUND)
    message(WARNING "WEST requires Python3 development libraries to be installed. Please install the appropriate package for your system (e.g., python3-dev or python3-devel) and re-run CMake.")
endif() 

# ==============================================================================
# STEP 4: compiling the modules and executables of WEST 
# ==============================================================================
message(STATUS "\n=======================================================")
message(STATUS ">>> STEP 4: Compiling WEST <<<")
message(STATUS "=======================================================\n")


#The order here is important, as the libraries are needed for the modules to compile correctly
add_subdirectory(Libraries/Forpy)
add_subdirectory(Libraries/Base64)
add_subdirectory(Libraries/Json)
add_subdirectory(Pytools)
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
add_subdirectory(Westpp)