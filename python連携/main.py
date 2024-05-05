import ctypes
import numpy as np
import sys

def call_fortran(n1, n2):
	dll = np.ctypeslib.load_library("main.dll", ".")
	dll.main.argtypes = [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32)]
	dll.main.restype = ctypes.c_void_p
	fn1 = ctypes.byref(ctypes.c_int32(n1))
	fn2 = ctypes.byref(ctypes.c_int32(n2))
	dll.TEST(fn1, fn2)

N1 = 61
N2 = 100
call_fortran(N1, N2)
print(N1)
