clear all;
close all;

x = linspace(-15,15,80);
[xx, yy] = meshgrid(x,x);

N = 20;
rng(2);
xc = [10*(2*rand(1,N)-1)];
yc = [10*(2*rand(1,N)-1)];

zz = zeros(size(xx));

sig2 = 8;

for n = [1:N]
	zz = zz + exp(-((xx-xc(n)).^2 + (yy-yc(n)).^2)/sig2);
end

surf(xx, yy, zz);

view(80,45);

print('-dpng', 'figure.png');
