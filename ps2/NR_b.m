dt=0.01; %timestep
Tfinal=10;
t=[0:dt:Tfinal]; %time vector
u0=[0; pi/4]; %initial condition
%BE method
U=zeros(2,length(t));
U(:,1)=u0;
g=9.8;L=1; %gravitational acc and pendulum length
for i=2:length(t)
    w = U(:,i-1); %start with an initial guess of the new value
    normR=100;
while normR > 1e-12 
    R = w - U(:,i-1) - dt*([-(g/L)*sin(w(2)); w(1)]); %Residual
    J = [0 -(g/L)*cos(w(2)); 1 0]; %Jacobian of the system 
    dw = (eye(2)-dt*J)\-R; % solve for dw 
    w = w + dw; 
    normR = norm(R) ;
end
U(:,i)=w;
end
%trapezoidal integration method
U_trp=zeros(2,length(t));
U_trp(:,1)=u0;
g=9.8;L=1; %gravitational acc and pendulum length
for i=2:length(t)
    w = U_trp(:,i-1); %start with an initial guess of the new value
    normR=100;
while normR > 1e-12 
    R = w - U_trp(:,i-1) - 0.5*dt*([-(g/L)*sin(w(2)); w(1)]+[-(g/L)*sin(U_trp(2,i-1)); U_trp(1,i-1)]); %Residual
    J = [0 -(g/L)*cos(w(2)); 1 0]; %Jacobian of the system 
    dw = (eye(2)-0.5*dt*J)\-R; % solve for dw 
    w = w + dw; 
    normR = norm(R) ;
end
U_trp(:,i)=w;
end

figure (1);
plot(t,U(2,:),'-',t,U_trp(2,:),'-')
legend('BE with NR','trapezoidal with NR');
title('Pendulum angular displacement calculated using BE/NR and trapezoidal/NR');
xlabel('time (seconds)');
ylabel('Angular displacement (radians)');
