dt=0.01; %timestep
u0=[0; pi/4]; %initial condition
%BE method
g=9.8;L=1; %gravitational acc and pendulum length
omega=sqrt(g/L);
Tfinal=2*pi*sqrt(L/g);
t=[0:dt:Tfinal]; %time vector
U=zeros(2,length(t),59);
U_analitic=zeros(2,length(t),59);

for indice=1:59 %initialize 3d matrix wit angular amplitude, knowing that initial ang speed is zero
    U(2,1,indice)=pi/(indice+1);
    U_analitic(2,1,indice)=pi/(indice+1);
end

for indice=1:59
    for i=2:length(t)
        w = U(:,i-1,indice); %start with an initial guess of the new value
        normR=100;
    while normR > 1e-12 
        R = w - U(:,i-1,indice) - dt*([-(g/L)*sin(w(2)); w(1)]); %Residual
        J = [0 -(g/L)*cos(w(2)); 1 0]; %Jacobian of the system 
        dw = (eye(2)-dt*J)\-R; % solve for dw 
        w = w + dw; 
        normR = norm(R) ;
    end
    U(:,i,indice)=w; %update matrix U at location i-1, indice with a vector
    U_analitic(:,i,indice)=[-U_analitic(2,1,indice)*omega*sin(omega*(i-1)*dt),U_analitic(2,1,indice)*cos(omega*(i-1)*dt)];
    end
end
relative_errors=zeros(1,59);
for indice=1:59 %look for the maximal absolute relative error
    %relative_errors(indice)=norm(U_analitic(:,length(t),indice)-U(:,length(t),indice))/norm(U_analitic(:,length(t),indice));
    relative_errors(indice)=norm(U_analitic(2,:,indice)-U(2,:,indice))/norm(U_analitic(2,:,indice));
end
disp(relative_errors)
figure (1);
plot(t,U(2,:,59),'-',t,U_analitic(2,:,59),'-')
legend('BE/NR','analitic');
title('Pendulum angular displacement');
xlabel('time (seconds)');
ylabel('Angular displacement (radians)');
