#include "fintrf.h"

module my_module
	implicit none
	private

	public :: parameter, timestwo

	type :: parameter
		real(kind=8) :: factor
	end type parameter

	contains

	subroutine timestwo(y, x, size, p)
		real(kind=8), intent(in)  :: x(*)
		real(kind=8), intent(out) :: y(*)
		type(parameter), intent(in) :: p
		integer n, size
		do n = 1, size
			y(n) = p%factor * x(n)
		end do

	end subroutine timestwo
end module my_module

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	use mpi
	use m_HACApK_base
	use m_HACApK_calc_entry_ij
	use my_module
	implicit none

	mwPointer plhs(*), prhs(*)
	integer nlhs, nrhs

	mwPointer mxGetDoubles
	mwPointer mxCreateDoubleMatrix

	integer mxIsNumeric
	mwPointer mxGetM, mxGetN

	mwPointer x_ptr, y_ptr
	mwPointer mrows, ncols
	mwSize size

	real(kind=8), allocatable :: x(:)
	real(kind=8), allocatable :: y(:)
	type(parameter) :: p

	mwPointer factor_ptr

	if(nrhs .ne. 2) then
		call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput', 'Two inputs required.')
	elseif(nlhs .ne. 1) then
		call mexErrMsgIdAndTxt ('MATLAB:timestwo:nOutput', 'One output required.')
	endif
	
	mrows = mxGetM(prhs(1))
	ncols = mxGetN(prhs(1))
	size = mrows*ncols

	allocate(x(size))
	allocate(y(size))

	x_ptr = mxGetDoubles(prhs(1))
	call mxCopyPtrToReal8(x_ptr, x, size)

	factor_ptr = mxGetDoubles(prhs(2))
	call mxCopyPtrToReal8(factor_ptr, p%factor, 1)

	call timestwo(y, x, size, p)

	plhs(1) = mxCreateDoubleMatrix(mrows, ncols, 0)
	y_ptr = mxGetDoubles(plhs(1))
	call mxCopyReal8ToPtr(y, y_ptr, size)
	
	deallocate(x)
	deallocate(y)
return
end

