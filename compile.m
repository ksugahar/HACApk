clear all;
close all;
setenv("ONEAPI_ROOT", "C:\Program Files (x86)\Intel\oneAPI")
mex -v -R2018a COMPFLAGS="$COMPFLAGS /free /QaxCORE-AVX2" -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"S:\Fortran\2024_08_31_mexの練習" main_mex.f90 m_HACApK_calc_entry_ij.obj HACApK_lib.obj

factor = 0.1;

L1 = 0.1000;
N1 = 61;

dLx1 = 2.0*L1/(N1-1);
dLy1 = 2.0*L1/(N1-1);

dx = [ -1.0, 1.0, 1.0,-1.0]*dLx1/2.0
dy = [ -1.0,-1.0, 1.0, 1.0]*dLy1/2.0
dz = [ -0.0, 0.0, 0.0, 0.0]

a = main_mex(reshape([1:63], 9, 7), factor, dx, dy, dz);
disp(a)

