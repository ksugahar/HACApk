clear all;
close all;
setenv("ONEAPI_ROOT", "C:\Program Files (x86)\Intel\oneAPI")
mex -v -R2018a COMPFLAGS="$COMPFLAGS /free /QaxCORE-AVX2" -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\mkl" main_mex.f90 HACApK_lib.obj m_HACApK_base.obj m_HACApK_calc_entry_ij.obj -L"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\lib" -L"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib" -llibmpi_ilp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
run_check



