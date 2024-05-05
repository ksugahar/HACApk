del /f main.dll
del *.o
del *.obj
del main.obj
del main.exp
del main.lib
ifx /dll /iface:cvf /nologo -o main.dll main.f90 HACApK_lib.f90 m_HACApK_base.f90 m_HACApK_calc_entry_ij.f90 mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib impi.lib libmpi_ilp64.lib
python main.py

