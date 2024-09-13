program main
use mpi
use m_HACApK_base
use m_HACApK_calc_entry_ij

integer*4 :: N1, N2
real*4 :: L1 = 0.1000d0, dLx1, dLy1, b = 0.059d0
real*4 :: L2 = 0.0015d0, dLx2, dLy2
integer*4 :: nx, ny ,nz
integer*4 :: i, j
type(st_HACApK_calc_entry) :: zbemv
type(st_HACApK_lcontrol) :: st_ctl
real*8, allocatable :: zab(:,:),zaa(:,:),param(:)
real*8 :: zeps, znrmmat, ACA_EPS
integer*4, allocatable :: lodl(:)
integer*4 :: nd,kmax,icomm,ierr
real*8,pointer, dimension(:) :: coil_x, coil_y, coil_z
real*8,pointer, dimension(:) :: eval_x, eval_y, eval_z
real*8,pointer, dimension(:) :: dx, dy, dz

!  各コイルの重心座標
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
	

!  測定点の座標

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
!!!Hz = loop_Biot_Savart(i,j)
Hz = HACApK_entry_ij(i,j,zbemv)
write(*,*) Hz

nd = max(n1,n2)
kmax=30

icomm = MPI_COMM_WORLD; ierr = 0
call MPI_Init ( ierr )
lrtrn=Hacapk_Init(nd,st_ctl,zbemv,icomm)
allocate(zaa(n1,kmax),zab(n2,kmax),lodl(nd))
zaa=0.0d0; zab=0.0d0

do il=1,nd;
	lodl(il)=il;
enddo

nstrtl=1; nstrtt=1

ACA_EPS=1.0e-15; znrmmat=1.0e-15

zeps=1.0e-12
kt = HACApK_acaplus(zaa,zab,st_ctl%param,n1,n2,nstrtl,nstrtt,lodl,zbemv,kmax,zeps,znrmmat,ACA_EPS)
print*,'Matrix rank; kt=',kt

! call Chk_rank_residual(zaa,zab,n1,n2,kt,kmax,zbemv)

end program main

subroutine Chk_rank_residual(zaa,zab,ndl,ndt,kt,kmax,zbemv)
use m_HACApK_calc_entry_ij
implicit real*8(a-h,o-z)
real*8 :: zab(1:ndt,1:kmax),zaa(1:ndl,1:kmax)
real*8,dimension(:),allocatable :: zresk
type(st_HACApK_calc_entry) :: zbemv
1000 format(5(a,i10))
2000 format(5(a,f12.2))

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
