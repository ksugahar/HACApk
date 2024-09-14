#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	use mpi
	use m_HACApK_base
	use m_HACApK_calc_entry_ij
	implicit none

	type(st_HACApK_calc_entry) :: zbemv
	type(st_HACApK_lcontrol) :: st_ctl
	real(kind=8), allocatable :: zab(:,:),zaa(:,:),param(:)
	real(kind=8) :: zeps, znrmmat, ACA_EPS
	integer(kind=4), allocatable :: lodl(:)
	integer(kind=4) :: nd, kmax, icomm, ierr

	mwPointer plhs(*), prhs(*)
	integer nlhs, nrhs

	mwPointer mxGetDoubles
	mwPointer mxCreateDoubleMatrix

	integer mxIsNumeric
	mwPointer mxGetM, mxGetN

	mwPointer n1, m1, n4, m4, n7, m7
	mwSize nm1, nm4, nm7

	real(kind=8) :: Hz
	integer(kind=4) :: i, j

	if(nrhs .ne. 9) then
		call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput', 'Nine inputs required.')
!	elseif(nlhs .ne. 1) then
!		call mexErrMsgIdAndTxt ('MATLAB:timestwo:nOutput', 'One output required.')
	endif
	
	n1 = mxGetN(prhs(1))
	m1 = mxGetM(prhs(1))
	nm1 = n1*m1
	allocate(zbemv%coil_x(nm1))
	allocate(zbemv%coil_y(nm1))
	allocate(zbemv%coil_z(nm1))
	call mxCopyPtrToReal8(mxGetDoubles(prhs(1)), zbemv%coil_x, nm1)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(2)), zbemv%coil_y, nm1)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(3)), zbemv%coil_z, nm1)

	n4 = mxGetN(prhs(4))
	m4 = mxGetM(prhs(4))
	nm4 = n4*m4
	allocate(zbemv%eval_x(nm4))
	allocate(zbemv%eval_y(nm4))
	allocate(zbemv%eval_z(nm4))
	call mxCopyPtrToReal8(mxGetDoubles(prhs(4)), zbemv%eval_x, nm4)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(5)), zbemv%eval_y, nm4)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(6)), zbemv%eval_z, nm4)

	n7 = mxGetN(prhs(7))
	m7 = mxGetM(prhs(7))
	nm7 = n7*m7
	allocate(zbemv%dx(nm7))
	allocate(zbemv%dy(nm7))
	allocate(zbemv%dz(nm7))
	call mxCopyPtrToReal8(mxGetDoubles(prhs(7)), zbemv%dx, nm7)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(8)), zbemv%dy, nm7)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(9)), zbemv%dz, nm7)

	i = 20
	j = 40
	Hz = HACApK_entry_ij(i, j, zbemv)

!	open(unit=1, file='result.txt', status='new', action='write')
!	write(1,*) Hz
!	close(unit=1, status='keep')

	plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
	call mxCopyReal8ToPtr(Hz, mxGetDoubles(plhs(1)), 1)

return
end

