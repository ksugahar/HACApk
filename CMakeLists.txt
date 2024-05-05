cmake_minimum_required(VERSION 3.10)
set(CMAKE_Fortran_COMPILER "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/ifx.exe")
project(Biot_Savart Fortran)

set(MKL_LIBRARIES "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/mkl_intel_lp64.lib"
									"C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/mkl_intel_thread.lib"
									"C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/mkl_core.lib"
									"C:/Program Files (x86)/Intel/oneAPI/compiler/latest/lib/libiomp5md.lib")
#include_directories("C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include")

set(MPI_LIBRARIES "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/lib/impi.lib"
									"C:/Program Files (x86)/Intel/oneAPI/mpi/latest/lib/libmpi_ilp64.lib")
include_directories("C:/Program Files (x86)/Intel/oneAPI/mpi/latest/include")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -shared")

add_library(main SHARED main.f90 HACApK_lib.f90 m_HACApK_base.f90 m_HACApK_calc_entry_ij.f90)
set_target_properties(main PROPERTIES OUTPUT_NAME "main")
target_link_directories(main PRIVATE "C:/Program Files (x86)/Intel/oneAPI/mkl")
target_link_libraries(main PRIVATE ${MKL_LIBRARIES} ${MPI_LIBRARIES})

