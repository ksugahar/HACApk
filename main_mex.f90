#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	use mpi
	use m_HACApK_base
	use m_HACApK_calc_entry_ij
	implicit none

	real(kind=8), allocatable :: C1(:,:), C2(:,:), SA(:,:), C4(:,:), Ca(:,:), Cb(:,:), dC(:,:)
	real(kind=8), allocatable :: C01(:,:), C02(:,:), SB(:,:)
	real(kind=8), dimension(3600,100) :: Hz
	real(kind=8), dimension(3600) :: I_HACApk
	real(kind=8), dimension(100) :: B_target
	integer(kind=4) :: nx, ny ,nz
	integer(kind=4) :: i, j, k
	integer(kind=4) :: kt, kc, kd
	integer(kind=4) :: lrtrn, nstrtt, nstrtl, nm01, nm04, il, it

	type(st_HACApK_calc_entry) :: zbemv
	type(st_HACApK_lcontrol) :: st_ctl
	real(kind=8), allocatable :: zab(:,:),zaa(:,:),param(:)
	type(st_svd) :: zac,zad
	real(kind=8) :: zeps,znrmmat,ACA_EPS,zeps1,zeps2
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

	if(nrhs .ne. 9) then
		call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput', 'Nine inputs required.')
	elseif(nlhs .ne. 1) then
		call mexErrMsgIdAndTxt ('MATLAB:timestwo:nOutput', 'One output required.')
	endif
	n1 = mxGetN(prhs(1))
	m1 = mxGetM(prhs(1))
	nm1 = n1*m1
	nm01 = n1 * m1
	allocate(zbemv%coil_x(nm1))
	allocate(zbemv%coil_y(nm1))
	allocate(zbemv%coil_z(nm1))
	call mxCopyPtrToReal8(mxGetDoubles(prhs(1)), zbemv%coil_x, nm1)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(2)), zbemv%coil_y, nm1)
	call mxCopyPtrToReal8(mxGetDoubles(prhs(3)), zbemv%coil_z, nm1)

	n4 = mxGetN(prhs(4))
	m4 = mxGetM(prhs(4))
	nm4 = n4*m4
	nm04 = n4*m4
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

	nd=max(nm01,n7)
	kmax=30
	icomm = MPI_COMM_WORLD; ierr = 0
	call MPI_Init ( ierr )
	lrtrn=Hacapk_Init(nd,st_ctl,zbemv,icomm)
	allocate(zaa(nm1,kmax),zab(nm4,kmax),lodl(nd))
	zaa=0.0d0; zab=0.0d0
	do il=1,nd; lodl(il)=il; enddo
	nstrtl=1; nstrtt=1
	ACA_EPS=1.0d-15; znrmmat=1.0d-15

	zeps=1.0e-12
	kt=HACApK_acaplus(zaa,zab,st_ctl%param,nm01,nm04,nstrtl,nstrtt,lodl,zbemv,kmax,zeps,znrmmat,ACA_EPS)

	! call Chk_rank_residual(zaa,zab,nm1,n7,kt,kmax,zbemv)
	zeps1=1.0d-5; zeps2=1.0d-4
	call TruncatedSVD(zeps1,zac,zaa,nm1,kt)
	call TruncatedSVD(zeps2,zad,zab,n7,kt)

	do i = 1, n7
		B_target(i) = zbemv%eval_x(i)
	end do

	kc = zac%k
	kd = zad%k
	allocate(C1(kt,nm1), C2(kt,nm1), SA(kt,nm1))

	do nx = 1, kc
		SA(nx,nx) = 1 / zac%zs(nx)
	end do

	do j = 1, size(SA,dim=2)
		do i = 1, size(zac%zv,dim=2)
			C1(i,j) = 0.0d0
			do k = 1, size(zac%zv,dim=1)
				C1(i,j) = C1(i,j) + zac%zv(k,i) * SA(k,j)
			end do
		end do
	end do

	do j = 1, size(zac%zu,dim=1)
		do i = 1, size(C1,dim=1)
			do k = 1, size(C1,dim=2)
				C2(i,j) = C2(i,j) + C1(i,k) * zac%zu(j,k)
			end do
		end do
	end do
	C2(nm1,kt) = 0.0d0

	allocate(C01(n7,kt), C02(n7,kt), SB(n7,kt), C4(n7,nm1))
	do nx = 1, kd
		SB(nx,nx) = 1 / zad%zs(nx)
	end do

	do j = 1, size(SB,dim=2)
		do i = 1, size(zad%zu,dim=1)
			C01(i,j) = 0.0d0
			do k = 1, size(zad%zu,dim=2)
				C01(i,j) = C01(i,j) + zad%zu(i,k) * SB(k,j)
			end do
		end do
	end do

	do j = 1, size(zad%zv,dim=2)
		do i = 1, size(C01,dim=1)
			do k = 1, size(C01,dim=2)
				C02(i,j) = C02(i,j) + C01(i,k) * zad%zv(k,j)
			end do
		end do
	end do
	C02(n7,kt) = 0.0d0

	do j = 1, size(C2,dim=2)
		do i = 1, size(C02,dim=1)
			C4(i,j) = 0.0d0
			do k = 1, size(C02,dim=2)
				C4(i,j) = C4(i,j) + C02(i,k) * C2(k,j)
			end do
		end do
	end do

	do i = 1, size(C4,dim=2)
		do j = 1, size(B_target,dim=1)
			I_HACApk(i) = I_HACApk(i) + B_target(j) * C4(j,i)
		end do
	end do

	plhs(1) = mxCreateDoubleMatrix(nm1, 1, 0)
	call mxCopyReal8ToPtr(I_HACApk, mxGetDoubles(plhs(1)), size(I_HACApk))

return
end


subroutine TruncatedSVD(zeps,zsvd,za,nn,k)
use m_HACApK_calc_entry_ij
	implicit real(kind=8)(a-h,o-z)
	real(kind=8) :: za(nn,k)
	type(st_svd) :: zsvd
	integer(kind=4),dimension(:), allocatable :: lrow_msk,lcol_msk
	real(kind=8),dimension(:),allocatable :: w,work
	real(kind=8),dimension(:,:),allocatable :: u,vt,waa,waa2

	1000 format(5(a,i10))
	2000 format(5(a,f12.2))
	3000 format(50(1pe11.2))

	lwork=10*nn; lda=nn; ldu=nn; ldvt=k
	allocate(w(nn),work(lwork),u(nn,nn),vt(nn,nn))
	call dgesvd ( 'S', 'S', nn, k, za, lda, w, u, ldu, vt, nn, work, lwork, info)

 !print*,'TruncatedSVD; zeps=',zeps
 !print*,'TruncatedSVD; e=',w(1:k+1)

	knw=k; zsvderr=0.0d0
	do il=2,k
		zsvderr=w(il)/w(1)
		if(zsvderr<zeps)then
			knw=il-1
			exit
		endif
	enddo
        !print*,'rank_e4=',knw,'/',nn
        !print*,'sv_max=',w(1)
        !print*,'sv_error=',w(knw+1)
 
 print*, 'TruncatedSVD; Truncated rank=',knw
 print*, 'TruncatedSVD; Truncated SVD error=',zsvderr
 zsvd%k=knw
 allocate(zsvd%zu(nn,nn),zsvd%zv(k,k),zsvd%zs(knw))
 zsvd%zs(:knw)=w(:knw)
 zsvd%zu(:nn,:nn)=u(:nn,:nn)
 zsvd%zv(:k,:k)=vt(:k,:k)
 
 return
 
     print*,'u'
     do il=1,nn
       w(:knw+1)=u(il,:knw+1)
       write(*,3000) w(:knw+1)
     enddo
     print*,'zsvd%zu'
     do il=1,nn
       w(:knw)=zsvd%zu(il,:knw)
       write(*,3000) w(:knw)
     enddo
     print*,'vt'
     do il=1,knw+1
       w(:k)=vt(il,:k)
       write(*,3000) w(:k)
     enddo
     print*,'zsvd%zv'
     do il=1,knw
       w(:k)=zsvd%zv(il,:k)
       write(*,3000) w(:k)
     enddo
 
 end subroutine

subroutine Chk_rank_residual(zaa,zab,ndl,ndt,kt,kmax,zbemv)
	use m_HACApK_calc_entry_ij
	implicit real(kind=8)(a-h,o-z)
	real(kind=8) :: zab(1:ndt,1:kmax),zaa(1:ndl,1:kmax)
	real(kind=8),dimension(:),allocatable :: zresk
	type(st_HACApK_calc_entry) :: zbemv
1000 format(5(a,i10))
2000 format(5(a,f12.2))
	integer(kind=4) :: il, it

 write(*,1000) 'ndl=',ndl,'; ndt=',ndt,'; kt=',kt
 allocate(zresk(0:kt)); zresk=0.0d0
 zanrm=0.0d0
 do il=1,ndl
   do it=1,ndt
     zz=HACApK_entry_ij(il,it,zbemv)
     zanrm=zanrm+zz*zz
     do ik=1,kt
       zz=zz-zaa(il,ik)*zab(it,ik)
       zresk(ik)=zresk(ik)+zz*zz
     enddo
   enddo
 enddo
 zresk(0)=zanrm
 do ik=0,kt
   zrr=dsqrt(zresk(ik)/zanrm)
   write(17,*) ik,zrr,log10(zrr)
 enddo
end subroutine
