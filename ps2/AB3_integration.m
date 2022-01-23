u0=1;
dt=0.01; %timestep
Tfinal=2; %finish time
t=[0:dt:Tfinal]; %time vector
U_AB=zeros(1,length(t));U_AB(1)=u0;
U_exact=zeros(1,length(t));U_exact(1)=u0;
for i=2:4
    U_AB(i)=U_AB(i-1)+dt*odefun(U_AB(i-1));
    U_exact(i)=exp(1)^(-5*i*dt);
end
for i=4:length(t)
    U_AB(i)=U_AB(i-1)+(dt/12)*(23*odefun(U_AB(i-1))-16*odefun(U_AB(i-2)) +5* odefun(U_AB(i-3)));
    U_exact(i)=exp(1)^(-5*i*dt);
end
figure(1);
plot(t,U_AB,'-',t,U_exact,'-')
legend('AB3','exact sol')
title('Solution comparison for dt=0.01');
xlabel('t (independent variable)');
ylabel('U');

function dudt=odefun(u)
dudt=-5*u
end