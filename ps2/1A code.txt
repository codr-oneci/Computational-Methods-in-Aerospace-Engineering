th = 0:pi/50:2*pi;
g=exp(1j*th);
z=zeros(1,length(g));
for i=1:length(g)
    c=g(i);
    z(i)=(12*c^3-12*c^2)/(23*c^2-16*c+5);
end
disp(g);
xunit = real(z);
yunit = imag(z);
h = plot(xunit, yunit);
title('Adams Bashforth with p=3 stability region ');
xlabel('Real \lambda\Delta t');
ylabel('Imag \lambda\Delta t');