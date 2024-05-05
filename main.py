import ctypes
import numpy as np
import sys

def call_fortran(n1, n2):
	dll = np.ctypeslib.load_library("main.dll", ".")
	dll.MAIN.argtypes = [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32)]
	dll.MAIN.restype = ctypes.c_void_p
	fn1 = ctypes.byref(ctypes.c_int32(n1))
	fn2 = ctypes.byref(ctypes.c_int32(n2))
	dll.MAIN(fn1, fn2)

N1 = 61
N2 = 100
call_fortran(N1, N2)
