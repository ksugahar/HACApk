clear all;
close all;

L1 = 0.1000;
L2 = 0.0015;
b = 0.059;
N1 = 61;
N2 = 100;

dLx1 = 2.0*L1/(N1-1);
dLy1 = 2.0*L1/(N1-1);
dLx2 = 2.0*L2/(N2-1);
dLy2 = 2.0*L2/(N2-1);

n = 0;
for ny = [1:N1]
	for nx = [1:N1]
		n = n + 1;
		coil_x(n) = -L1 + (nx-1)*dLx1;
		coil_y(n) = -L1 + (ny-1)*dLy1;
		coil_z(n) = b/2;
	end
end
n = 0;
for ny = [1:N2]
	for nx = [1:N2]
		n = n + 1;
		eval_x(n) = -L2 + (nx-1)*dLx2;
		eval_y(n) = -L2 + (ny-1)*dLy2;
		eval_z(n) = 0;
	end
end

dx = [ -1.0, 1.0, 1.0,-1.0]*dLx1/2.0;
dy = [ -1.0,-1.0, 1.0, 1.0]*dLy1/2.0;
dz = [ -0.0, 0.0, 0.0, 0.0];

a = main_mex(coil_x, coil_y, coil_z, eval_x, eval_y, eval_z, dx, dy, dz);

format long
disp(a)
