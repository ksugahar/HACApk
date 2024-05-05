del /f main.dll
del *.o
del *.obj
ifx -c HACApK_lib.f90 -o HACApK_lib.obj
ifx -c m_HACApK_base.f90 -o m_HACApK_base.obj
ifx -c m_HACApK_calc_entry_ij.f90 -o m_HACApK_calc_entry_ij.obj
ifx -c main.f90 -o main.obj
rem ifx /Qopenmp /QaxCORE-AVX2 main.obj m_HACApK_calc_entry_ij.obj m_HACApK_base.f90 HACApK_lib.f90 mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib /dll /out:main.dll 
ifx /Qopenmp /QaxCORE-AVX2 main.obj m_HACApK_calc_entry_ij.obj m_HACApK_base.f90 HACApK_lib.f90 mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib /dll -static
python main.py
del main.exp
del main.lib

