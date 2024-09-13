del *.obj
del *.mod
del main.exe
ifx /Qopenmp /QaxCORE-AVX2 -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"S:\Fortran\2024_08_31_mex‚Ì—ûK" -c HACApK_lib.f90               mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib
ifx /Qopenmp /QaxCORE-AVX2 -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"S:\Fortran\2024_08_31_mex‚Ì—ûK" -c m_HACApK_calc_entry_ij.f90   mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib
ifx /Qopenmp /QaxCORE-AVX2 -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"S:\Fortran\2024_08_31_mex‚Ì—ûK" -c m_HACApK_base.f90            mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib
ifx /Qopenmp /QaxCORE-AVX2 -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include\mpi" -I"S:\Fortran\2024_08_31_mex‚Ì—ûK" /exe:main main.f90 m_HACApK_base.obj m_HACApK_calc_entry_ij.obj HACApK_lib.obj mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib
main.exe
