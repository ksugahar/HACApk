clear all;
close all;
setenv("ONEAPI_ROOT", "C:\Program Files (x86)\Intel\oneAPI")
mex -v -R2018a COMPFLAGS="$COMPFLAGS /free /QaxCORE-AVX2" -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"S:\Fortran\2024_08_31_mexの練習" main_mex.f90 m_HACApK_calc_entry_ij.obj HACApK_lib.obj

run_check

