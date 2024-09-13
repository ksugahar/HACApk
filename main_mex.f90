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

	integer(kind=4) :: N1, N2
	real(kind=8) :: L1 = 0.1000d0, dLx1, dLy1, b = 0.059d0
	real(kind=8) :: L2 = 0.0015d0, dLx2, dLy2
	integer(kind=4) :: nx, ny ,nz
	integer(kind=4) :: i, j
	type(st_HACApK_calc_entry) :: zbemv
	type(st_HACApK_lcontrol) :: st_ctl
	real(kind=8), allocatable :: zab(:,:),zaa(:,:),param(:)
	real(kind=8) :: zeps, znrmmat, ACA_EPS
	integer(kind=4), allocatable :: lodl(:)
	integer(kind=4) :: nd,kmax,icomm,ierr
	real(kind=8), pointer, dimension(:) :: coil_x, coil_y, coil_z
	real(kind=8), pointer, dimension(:) :: eval_x, eval_y, eval_z
	real(kind=8), pointer, dimension(:) :: dx, dy, dz
	real(kind=8) :: Hz

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

	N1 = 61
	allocate(zbemv%coil_x(N1*N1))
	allocate(zbemv%coil_y(N1*N1))
	allocate(zbemv%coil_z(N1*N1))
	coil_x=>zbemv%coil_x(:); coil_y=>zbemv%coil_y; coil_z=>zbemv%coil_z

	dLx1 = 2.0d0*L1/dble(N1-1)
	dLy1 = 2.0d0*L1/dble(N1-1)
	i = 0
	do ny = 1, N1
		do nx = 1, N1
			i = i + 1
			coil_x(i) = -L1 + dble(nx-1)*dLx1
			coil_y(i) = -L1 + dble(ny-1)*dLy1
			coil_z(i) = b/2
		end do
	end do

	N2 = 100
	allocate(zbemv%eval_x(N2*N2))
	allocate(zbemv%eval_y(N2*N2))
	allocate(zbemv%eval_z(N2*N2))
	eval_x=>zbemv%eval_x; eval_y=>zbemv%eval_y; eval_z=>zbemv%eval_z
	dLx2 = 2.0d0*L2/dble(N2-1)
	dLy2 = 2.0d0*L2/dble(N2-1)
	j = 0
	do ny = 1, N2
		do nx = 1, N2
			j = j + 1
			eval_x(j) = -L2 + dble(nx-1)*dLx2
			eval_y(j) = -L2 + dble(ny-1)*dLy2
			eval_z(j) = 0
		end do
	end do

	allocate(zbemv%dx(4))
	allocate(zbemv%dy(4))
	allocate(zbemv%dz(4))
	dx=>zbemv%dx(:); dy=>zbemv%dy; dz=>zbemv%dz

	dx = (/ -1.d0, 1.d0, 1.d0,-1.d0/)*dLx1/2.0
	dy = (/ -1.d0,-1.d0, 1.d0, 1.d0/)*dLy1/2.0
	dz = (/ -0.d0, 0.d0, 0.d0, 0.d0/)

	i = 20
	j = 40
	Hz = HACApK_entry_ij(i, j, zbemv)
	write(1,*) Hz

	plhs(1) = mxCreateDoubleMatrix(mrows, ncols, 0)
	y_ptr = mxGetDoubles(plhs(1))
	call mxCopyReal8ToPtr(y, y_ptr, size)
	
	deallocate(x)
	deallocate(y)
return
end

