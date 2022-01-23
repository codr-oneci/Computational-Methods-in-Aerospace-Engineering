% Specify x range and number of points 
x0 = -10; x1 = 10; 
Nx = 1000;
% Specify y range and number of points 
y0 = -10; y1 = 10; 
Ny = 1000;
% Construct mesh
xv = linspace(x0,x1,Nx);
yv = linspace(y0,y1,Ny);
[x,y] = meshgrid(xv,yv);
% Calculate z 
z = x + i*y;
% 2nd order Runge-Kutta growth factor 
%g = 1 + z + 0.5*z.^2;
a=12;b=-12-23*z;c=16*z;d=-5*z;
p=-b/(3*a);q=p^3+(b*c-3*a*d)/(6*a^2); r=c/(3*a);
g=(q+(((q^2)+(r-p^2)^3)^(1/2)))^(1/3)+(q-(((q^2)+(r-p^2)^3)^(1/2)))^(1/3)+p
% Calculate magnitude of g 
gmag = abs(g);
% Plot contours of gmag 
contour(x,y,gmag,[1 1],'k-');
axis([x0,x1,y0,y1]); 
axis('square'); 
xlabel('Real \lambda\Delta t');
ylabel('Imag \lambda\Delta t');
grid on;