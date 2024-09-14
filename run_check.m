clear all;
close all;
format long

L1 = 0.1000;
N1 = 61;
b = 0.059;

L2 = 0.0015;
N2 = 100;

[yy, xx] = meshgrid(linspace(-L1, L1, N1));
coil_x = xx(:);
coil_y = yy(:);
coil_z = b/2*ones(size(coil_x));

[yy, xx] = meshgrid(linspace(-L2, L2, N2));
eval_x = xx(:);
eval_y = yy(:);
eval_z = 0*ones(size(eval_x));

dLx1 = 2.0*L1/(N1-1);
dLy1 = 2.0*L1/(N1-1);

dx = [ -1.0, 1.0, 1.0,-1.0]*dLx1/2.0;
dy = [ -1.0,-1.0, 1.0, 1.0]*dLy1/2.0;
dz = [ -0.0, 0.0, 0.0, 0.0];

a = main_mex(coil_x, coil_y, coil_z, eval_x, eval_y, eval_z, dx, dy, dz);

disp(a)

