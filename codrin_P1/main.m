% Author: Codrin Oneci
%set constants that will be used in simulations (function setconst)

%SOLVE QUESTION 1
%Define function that gives the derivative of the state vector with regard
%to the independent variable (time t). See def at the end of script

% initialize state vector
u0=[22;0;0;5];
Opt= odeset('Events', @ground_intersection,'RelTol',1e-15,'AbsTol',1e-17,'NormControl','on'); %set options for ode45
[t,u] = ode45(@odefun,[0,100],u0,Opt); %solve for u
%%plot(t,u(:,3),'-',t,u(:,4),'-')
%plot glider's location on X and Y on the same figure
figure(1);
plot(u(:,3),u(:,4),'-')
title('Position Solution of fourth order glider system with ODE45');
xlabel('X location (m)');
ylabel('Y location (m)');
legend('ode45');
%plot glider's speed and angle of attack on the same figure
figure(2);
plot(u(:,2),u(:,1),'-')
title('Speed-Angle of Attack Solution of fourth order glider system with ODE45');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');
legend('ode45');
%g=sprintf('%d ', t);
%fprintf('Answer: %s\n', g)
disp(t(length(t))); %see stopping time here

%SOLVE QUESTION 2
%FORWARD EULER
dt=0.01; %timestep
Tfinal=t(length(t)); %time length of experiment
t=[0:dt:Tfinal]; %time vector
X=zeros(1,length(t));Y=zeros(1,length(t));V=zeros(1,length(t));AOA=zeros(1,length(t));
X(1)=u0(3);Y(1)=u0(4);AOA(1)=u0(2);V(1)=u0(1);
for i=2:length(t)
    state_derivative=odefun(i,[V(i-1);AOA(i-1);X(i-1);Y(i-1)]);
    V(i)=V(i-1)+dt*state_derivative(1);
    AOA(i)=AOA(i-1)+dt*state_derivative(2);
    X(i)=X(i-1)+dt*state_derivative(3);
    Y(i)=Y(i-1)+dt*state_derivative(4);
end
%plot glider's location on X and Y on the same figure
figure(3);
%plot(X,Y,'-');legend('FE');
%plot(u(:,3),u(:,4),'-');legend('ode45');
plot(X,Y,'-',u(:,3),u(:,4),'-')
legend('FE','ode45')
title('Position Solution of fourth order glider system with ODE45 and FE');
xlabel('X location (m)');
ylabel('Y location (m)');
%%legend('x','y')
%plot glider's speed and angle of attack on the same figure
figure(4);
%plot(AOA,V,'-');legend('FE');
%plot(u(:,2),u(:,1),'-');legend('ode45');
plot(AOA,V,'-',u(:,2),u(:,1),'-')
legend('FE','ode45');
title('Speed-Angle of Attack Solution of fourth order glider system with ODE45 and FE');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');


%RUNGE-KUTTA
dt=1.2; %timestep
Tfinal=t(length(t)); %time length of experiment
t=[0:dt:Tfinal]; %time vector
X=zeros(1,length(t));Y=zeros(1,length(t));V=zeros(1,length(t));AOA=zeros(1,length(t));
X(1)=u0(3);Y(1)=u0(4);AOA(1)=u0(2);V(1)=u0(1);
for i=2:length(t)
    vn=[V(i-1);AOA(i-1);X(i-1);Y(i-1)];
    a=dt*odefun(i*dt,vn);
    b=dt*odefun(i*dt+dt/2,vn+a/2);
    c=dt*odefun(i*dt+dt/2,vn+b/2);
    d=dt*odefun(i*dt+dt,vn+c);
    vn_new=vn+(a+2*b+2*c+d)/6;
    V(i)=vn_new(1);
    AOA(i)=vn_new(2);
    X(i)=vn_new(3);
    Y(i)=vn_new(4);
end
%plot glider's location on X and Y on the same figure
figure(5);
%plot(X,Y,'-');legend('FE');
%plot(u(:,3),u(:,4),'-');legend('ode45');
plot(X,Y,'-',u(:,3),u(:,4),'-')
legend('Runge-Kutta','ode45')
title('Position Solution of fourth order glider system with ODE45 and Runge-Kutta');
xlabel('X location (m)');
ylabel('Y location (m)');
%%legend('x','y')
%plot glider's speed and angle of attack on the same figure
figure(6);
%plot(AOA,V,'-');legend('FE');
%plot(u(:,2),u(:,1),'-');legend('ode45');
plot(AOA,V,'-',u(:,2),u(:,1),'-')
legend('Runge-Kutta','ode45');
title('Speed-Angle of Attack Solution of fourth order glider system with ODE45 and Runge-Kutta');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');



%SOLVE QUESTION 3 order of accuracy
%FE ACCURACY
u0=[22;0;0;5];
Opt= odeset('Events', @ground_intersection,'RelTol',1e-15,'AbsTol',1e-17,'NormControl','on'); %set options for ode45
[t_exact,u_exact] = ode45(@odefun,[0,10],u0,Opt); %solve for u
u_exact=u_exact.'.'.';
t_exact=t_exact.'.'.';
FE_errors=[0;0;0;0];
timesteps=[0.05;0.01;0.005;0.001];
%timesteps=[1;0.1;0.01;0.001;0.0001;0.00001];
for ts=1:length(timesteps)
    max_err=0;
    dt=timesteps(ts);
    Tfinal=10; %time length of experiment
    t=[0:dt:Tfinal]; %time vector
    V=zeros(4,length(t));
    V(:,1)=u0;
    for i=2:length(t) %FE iteration
        V(:,i)=V(:,i-1)+dt*odefun(i,V(:,i-1));
    end
    %disp(V);
    %create solution array of the same length as V
    %for i=1:length(t)
     %   ct=(i-1)*dt
     for i=2:length(t)
        [d,it]=min(abs(t_exact-i*dt)); %find element in t_exact closest to i*dt 
        err=norm(V(:,i)-u_exact(:,it));
        if err>max_err
         max_err=err;
     end
     end
    %end
    %disp(u_exact(1,:)-V(:,1))
     
     
     %disp (V(:,i)-u(it));
     %err=norm(V(:,i)-u_exact(it));
     %err=max(abs(V(:,i)-u_exact(it)));
     
    FE_errors(ts)=max_err; %maximal error measured in 10s
end

%RK ACCURACY
RK_errors=[0;0;0;0];
for ts=1:length(timesteps)
    max_err=0;
    dt=timesteps(ts);
    Tfinal=10; %time length of experiment
    t=[0:dt:Tfinal]; %time vector
    V=zeros(4,length(t));
    V(:,1)=u0;
    for i=2:length(t) %RK iteration
        vn=V(:,i-1);
    a=dt*odefun(i*dt,vn);
    b=dt*odefun(i*dt+dt/2,vn+a/2);
    c=dt*odefun(i*dt+dt/2,vn+b/2);
    d=dt*odefun(i*dt+dt,vn+c);
    V(:,i)=vn+(a+2*b+2*c+d)/6;
    end
    %disp(V);
    %create solution array of the same length as V
    %for i=1:length(t)
     %   ct=(i-1)*dt
     for i=2:length(t)
        [d,it]=min(abs(t_exact-i*dt)); %find element in t_exact closest to i*dt 
        err=norm(V(:,i)-u_exact(:,it));
        if err>max_err
         max_err=err;
     end
     end
    %end
    %disp(u_exact(1,:)-V(:,1))
     
     
     %disp (V(:,i)-u(it));
     %err=norm(V(:,i)-u_exact(it));
     %err=max(abs(V(:,i)-u_exact(it)));
     
    RK_errors(ts)=max_err; %maximal error measured in 10s
end
mem_V=V;

%disp(FE_errors);
figure (7)
plot(log(timesteps),log(FE_errors),'-',log(timesteps),log(RK_errors),'-')
legend('Forward Euler','Runge-Kutta 4');
title('Log Global Error - Log Timestep graph for order of accuracy determination');
xlabel('log timestep');
ylabel('log global error');


%PART 4
ua=[29;0;0;10];
ub=[23.1;0;0;10];
uc=[12.0223;-0.0831;0;10];
ud=[6;0;0;10];
dt=0.1;

%A
u0=ua;
[t,u] = ode45(@odefun,[0,100],u0,Opt); %use this to have a comparison; useful also for obtaining the specific end time
Tfinal=t(length(t)); %time length of experiment
t=[0:dt:Tfinal]; %time vector
X=zeros(1,length(t));Y=zeros(1,length(t));V=zeros(1,length(t));AOA=zeros(1,length(t));
X(1)=u0(3);Y(1)=u0(4);AOA(1)=u0(2);V(1)=u0(1);
for i=2:length(t)
    vn=[V(i-1);AOA(i-1);X(i-1);Y(i-1)];
    a=dt*odefun(i*dt,vn);
    b=dt*odefun(i*dt+dt/2,vn+a/2);
    c=dt*odefun(i*dt+dt/2,vn+b/2);
    d=dt*odefun(i*dt+dt,vn+c);
    vn_new=vn+(a+2*b+2*c+d)/6;
    V(i)=vn_new(1);
    AOA(i)=vn_new(2);
    X(i)=vn_new(3);
    Y(i)=vn_new(4);
end
%plot glider's location on X and Y on the same figure
figure(8);
%plot(X,Y,'-');legend('FE');
%plot(u(:,3),u(:,4),'-');legend('ode45');
plot(X,Y,'-',u(:,3),u(:,4),'-')
legend('Runge-Kutta','ode45')
title('Position Solution of fourth order glider system situation A');
xlabel('X location (m)');
ylabel('Y location (m)');
%%legend('x','y')
%plot glider's speed and angle of attack on the same figure
figure(9);
%plot(AOA,V,'-');legend('FE');
%plot(u(:,2),u(:,1),'-');legend('ode45');
plot(AOA,V,'-',u(:,2),u(:,1),'-')
legend('Runge-Kutta','ode45');
title('Speed-Angle of Attack Solution of fourth order glider system situation A');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');

%B
u0=ub;
[t,u] = ode45(@odefun,[0,100],u0,Opt); %use this to have a comparison; useful also for obtaining the specific end time
Tfinal=t(length(t)); %time length of experiment
t=[0:dt:Tfinal]; %time vector
X=zeros(1,length(t));Y=zeros(1,length(t));V=zeros(1,length(t));AOA=zeros(1,length(t));
X(1)=u0(3);Y(1)=u0(4);AOA(1)=u0(2);V(1)=u0(1);
for i=2:length(t)
    vn=[V(i-1);AOA(i-1);X(i-1);Y(i-1)];
    a=dt*odefun(i*dt,vn);
    b=dt*odefun(i*dt+dt/2,vn+a/2);
    c=dt*odefun(i*dt+dt/2,vn+b/2);
    d=dt*odefun(i*dt+dt,vn+c);
    vn_new=vn+(a+2*b+2*c+d)/6;
    V(i)=vn_new(1);
    AOA(i)=vn_new(2);
    X(i)=vn_new(3);
    Y(i)=vn_new(4);
end
%plot glider's location on X and Y on the same figure
figure(10);
%plot(X,Y,'-');legend('FE');
%plot(u(:,3),u(:,4),'-');legend('ode45');
plot(X,Y,'-',u(:,3),u(:,4),'-')
legend('Runge-Kutta','ode45')
title('Position Solution of fourth order glider system situation B');
xlabel('X location (m)');
ylabel('Y location (m)');
%%legend('x','y')
%plot glider's speed and angle of attack on the same figure
figure(11);
%plot(AOA,V,'-');legend('FE');
%plot(u(:,2),u(:,1),'-');legend('ode45');
plot(AOA,V,'-',u(:,2),u(:,1),'-')
legend('Runge-Kutta','ode45');
title('Speed-Angle of Attack Solution of fourth order glider system situation B');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');

%C
u0=uc;
[t,u] = ode45(@odefun,[0,100],u0,Opt); %use this to have a comparison; useful also for obtaining the specific end time
Tfinal=t(length(t)); %time length of experiment
t=[0:dt:Tfinal]; %time vector
X=zeros(1,length(t));Y=zeros(1,length(t));V=zeros(1,length(t));AOA=zeros(1,length(t));
X(1)=u0(3);Y(1)=u0(4);AOA(1)=u0(2);V(1)=u0(1);
for i=2:length(t)
    vn=[V(i-1);AOA(i-1);X(i-1);Y(i-1)];
    a=dt*odefun(i*dt,vn);
    b=dt*odefun(i*dt+dt/2,vn+a/2);
    c=dt*odefun(i*dt+dt/2,vn+b/2);
    d=dt*odefun(i*dt+dt,vn+c);
    vn_new=vn+(a+2*b+2*c+d)/6;
    V(i)=vn_new(1);
    AOA(i)=vn_new(2);
    X(i)=vn_new(3);
    Y(i)=vn_new(4);
end
%plot glider's location on X and Y on the same figure
figure(12);
%plot(X,Y,'-');legend('FE');
%plot(u(:,3),u(:,4),'-');legend('ode45');
plot(X,Y,'-',u(:,3),u(:,4),'-')
legend('Runge-Kutta','ode45')
title('Position Solution of fourth order glider system situation C');
xlabel('X location (m)');
ylabel('Y location (m)');
%%legend('x','y')
%plot glider's speed and angle of attack on the same figure
figure(13);
%plot(AOA,V,'-');legend('FE');
%plot(u(:,2),u(:,1),'-');legend('ode45');
plot(AOA,V,'-',u(:,2),u(:,1),'-')
legend('Runge-Kutta','ode45');
title('Speed-Angle of Attack Solution of fourth order glider system situation C');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');


%D
u0=ud;
[t,u] = ode45(@odefun,[0,100],ud,Opt); %use this to have a comparison; useful also for obtaining the specific end time
Tfinal=t(length(t)); %time length of experiment
t=[0:dt:Tfinal]; %time vector
X=zeros(1,length(t));Y=zeros(1,length(t));V=zeros(1,length(t));AOA=zeros(1,length(t));
X(1)=u0(3);Y(1)=u0(4);AOA(1)=u0(2);V(1)=u0(1);
for i=2:length(t)
    vn=[V(i-1);AOA(i-1);X(i-1);Y(i-1)];
    a=dt*odefun(i*dt,vn);
    b=dt*odefun(i*dt+dt/2,vn+a/2);
    c=dt*odefun(i*dt+dt/2,vn+b/2);
    d=dt*odefun(i*dt+dt,vn+c);
    vn_new=vn+(a+2*b+2*c+d)/6;
    V(i)=vn_new(1);
    AOA(i)=vn_new(2);
    X(i)=vn_new(3);
    Y(i)=vn_new(4);
end
%plot glider's location on X and Y on the same figure
figure(14);
%plot(X,Y,'-');legend('FE');
%plot(u(:,3),u(:,4),'-');legend('ode45');
plot(X,Y,'-',u(:,3),u(:,4),'-')
legend('Runge-Kutta','ode45')
title('Position Solution of fourth order glider system situation D');
xlabel('X location (m)');
ylabel('Y location (m)');
%%legend('x','y')
%plot glider's speed and angle of attack on the same figure
figure(15);
%plot(AOA,V,'-');legend('FE');
%plot(u(:,2),u(:,1),'-');legend('ode45');
plot(AOA,V,'-',u(:,2),u(:,1),'-')
legend('Runge-Kutta','ode45');
title('Speed-Angle of Attack Solution of fourth order glider system situation D');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');


%PART 6

c=setconst();
[x,y] = meshgrid(-1*pi:3*pi,1:30);
%disp(x);disp(size(x));
dim=size(x);
u = zeros(size(x));
v = zeros(size(x));
for th=1:dim(1)
    for speed=1:dim(2)
        u(th,speed)=c.RL*speed-c.g*cos(th)/(speed);
        v(th,speed)=-c.g*sin(th)-c.RD*speed^2;
    end
end


figure(16)
quiver(x,y,u,v)
title('2D vector respresentation on state change in AOA-V phase space');
xlabel('Angle of attack (rad)');
ylabel('Speed (m/s)');








%Define function that gives the derivative of the state vector with regard
%to the independent variable (time t)
function dudt = odefun(t,u)
Data=setconst();
dudt = [-Data.g*sin(u(2))-Data.RD*u(1)^2; Data.RL*u(1)-Data.g*cos(u(2))/u(1); u(1)*cos(u(2)); u(1)*sin(u(2)) ];
end


function [position, isterminal, direction] = ground_intersection(t, u)
 position = u(4); % event triggers when position = 0
 isterminal = 1; % halt integration
 direction = -1; % trigger when event function is decreasing
end

function [Data]=setconst()
%set constants that will be used in simulations
Data.g=9.81; %m/s^2
Data.air_density=1.22; %kg/m^3
Data.m=0.65; %kg
Data.S=0.06; %m^2
Data.CD=0.10; % drag coefficient (non-dimensional)
Data.CL=1.20; %lift coefficient (non-dimensional)
Data.RD=Data.air_density*Data.CD*Data.S/(2*Data.m); %physical constant used in the ODE and defined in problem's statement
Data.RL=Data.air_density*Data.CL*Data.S/(2*Data.m); %physical constant used in the ODE and defined in problem's statement
end

