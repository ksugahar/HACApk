module m_HACApK_calc_entry_ij

!*** type :: st_HACApK_calc_entry
	type :: st_HACApK_calc_entry
	real*8,pointer :: ao(:)
	integer :: nd,lp61
	
	real*8,pointer,dimension(:) :: coil_x, coil_y, coil_z
	real*8,pointer,dimension(:) :: eval_x, eval_y, eval_z
	real*8,pointer,dimension(:) :: dx, dy, dz
	
	end type st_HACApK_calc_entry
	
!*** type :: st_biot_physvar
	type :: st_biot_physvar
	real*8,pointer,dimension(:) :: coil_x, coil_y, coil_z
	real*8,pointer,dimension(:) :: eval_x, eval_y, eval_z
	real*8,pointer,dimension(:) :: dx, dy, dz
	end type st_biot_physvar 

public :: HACApK_entry_ij

contains
!***HACApK_entry_ij
	real*8 function HACApK_entry_ij(i,j,zbemv)
	integer*4 :: i, j
	type(st_HACApK_calc_entry) :: zbemv
	HACApK_entry_ij = loop_Biot_Savart(i,j,zbemv)
end function HACApK_entry_ij

	real*8 function loop_Biot_Savart(i, j,zbemv)
	real*8 :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4
	real*8, dimension(3) :: H12, H23, H34, H41, H
	integer*4 :: i, j
	type(st_HACApK_calc_entry) :: zbemv
	real*8,pointer,dimension(:) :: coil_x, coil_y, coil_z
	real*8,pointer,dimension(:) :: eval_x, eval_y, eval_z
	real*8,pointer,dimension(:) :: dx, dy, dz
	coil_x=>zbemv%coil_x(:); coil_y=>zbemv%coil_y; coil_z=>zbemv%coil_z
	eval_x=>zbemv%eval_x; eval_y=>zbemv%eval_y; eval_z=>zbemv%eval_z
	dx=>zbemv%dx(:); dy=>zbemv%dy; dz=>zbemv%dz
	X1 = coil_x(i) + dx(1)
	X2 = coil_x(i) + dx(2)
	X3 = coil_x(i) + dx(3)
	X4 = coil_x(i) + dx(4)
	Y1 = coil_y(i) + dy(1)
	Y2 = coil_y(i) + dy(2)
	Y3 = coil_y(i) + dy(3)
	Y4 = coil_y(i) + dy(4)
	Z1 = coil_z(i) + dz(1)
	Z2 = coil_z(i) + dz(2)
	Z3 = coil_z(i) + dz(3)
	Z4 = coil_z(i) + dz(4)

	!write(*,*) eval_x(j),eval_y(j),eval_z(j), X1, Y1, Z1, X2, Y2, Z2

	H12 = Biot_Savart(eval_x(j), eval_y(j), eval_z(j), X1, Y1, Z1, X2, Y2, Z2)
	H23 = Biot_Savart(eval_x(j), eval_y(j), eval_z(j), X2, Y2, Z2, X3, Y3, Z3)
	H34 = Biot_Savart(eval_x(j), eval_y(j), eval_z(j), X3, Y3, Z3, X4, Y4, Z4)
	H41 = Biot_Savart(eval_x(j), eval_y(j), eval_z(j), X4, Y4, Z4, X1, Y1, Z1)
	H = H12 + H23 + H34 + H41

	loop_Biot_Savart = H(3)
	end function

	function Biot_Savart(Ox, Oy, Oz, X1, Y1, Z1, X2, Y2, Z2)
		real*8 Biot_Savart(3)
		real*8 :: X1, Y1, Z1, X2, Y2, Z2
		real*8 :: Ox, Oy, Oz
		real*8 :: X21, Y21, Z21, R21
		real*8 :: OxX1, OyY1, OzZ1, OR1,OxX2, OyY2, OzZ2, OR2
		real*8 :: O121, O221, RcosO12, RcosO21, Lx1, Ly1, Lz1, L1
	real*8, parameter :: pi = 3.14159265358979323846d0
	real*8, parameter :: eps = 1d-15
	

		X21 = X2 - X1
		Y21 = Y2 - Y1
		Z21 = Z2 - Z1
		R21 = X21*X21 + Y21*Y21 + Z21*Z21

		OxX1 = Ox - X1
		OyY1 = Oy - Y1
		OzZ1 = Oz - Z1
		OR1 = SQRT(OxX1*OxX1 + OyY1*OyY1 + OzZ1*OzZ1)

		OxX2 = Ox - X2
		OyY2 = Oy - Y2
		OzZ2 = Oz - Z2
		OR2 = SQRT(OxX2*OxX2 + OyY2*OyY2 + OzZ2*OzZ2)

		O121 = OxX1*X21 + OyY1*Y21 + OzZ1*Z21
		O221 = OxX2*X21 + OyY2*Y21 + OzZ2*Z21

		RcosO12 = O121 / OR1
		RcosO21 = -O221 / OR2

		Lx1 = OxX1 - O121*X21 / R21
		Ly1 = OyY1 - O121*Y21 / R21
		Lz1 = OzZ1 - O121*Z21 / R21
		L1 = Lx1*Lx1 + Ly1*Ly1 + Lz1*Lz1 + eps

		Biot_Savart(1) = (Y21*Lz1 - Z21*Ly1) / (4.0d0*pi*L1*R21) * (RcosO12 + RcosO21)
		Biot_Savart(2) = (Z21*Lx1 - X21*Lz1) / (4.0d0*pi*L1*R21) * (RcosO12 + RcosO21)
		Biot_Savart(3) = (X21*Ly1 - Y21*Lx1) / (4.0d0*pi*L1*R21) * (RcosO12 + RcosO21)

	end function Biot_Savart

endmodule m_HACApK_calc_entry_ij
