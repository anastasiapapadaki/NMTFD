% The explicit Euler Method
clc
clear all
close all
% Parameters
T= 0.16;
dt=0.0001;
t= T/dt ;
dx=1/40;
dy=1/40;
Lx=1.0;
Ly=1.0;
Nx=(Lx/dx)+1;
Ny=(Ly/dy)+1;

% Space discretization
x = 0:dx:1;
y = 0:dy:1;

%Initial Condition
w=zeros(Nx,Ny,t); % we start with the initial temperature of zero everywhere

% Boundary conditions

for l=1:1:Ny % w(0, y, t) = 1-y3
  w(l,1,:) = 1-x(l)^3;
end

for r=1:1:Nx % w(1, y, t) = 1-sin(p/2 *y)
  w(r,Nx,:) = 1-sin(pi/2 * x(r));
end

w(1,:,:) = 1; % w(x, 0, t) = 1
w(Ny,:,:) = 0; % w(x, 1, t) = 0

% Discretization
for k=1:1:t; % time discretization
  k
  for i=2:1:Nx-1
    for j=2:1:Ny-1
      w(i,j,k+1) = w(i,j,k)+ dt *(((w(i-1,j,k)-2*w(i,j,k)+w(i+1,j,k))/dx^2)+
      ((w(i,j-1,k)-2*w(i,j,k)+w(i,j+1,k))/dy^2));
    endfor
  endfor
end

% Plotting the Contour_Plot
figure 1
w(:,:,1) % just to see the 1st iteration
contourf(x,y,w(:,:,t), 10);
colorbar;
xlabel('x (dimensionless)')
ylabel('y (dimensionless)')
title('Temperature at t=0.16s')
saveas(1,'Contour_Plot.png')

% Plotting the time evolution of the temperature at x=y=0.4 and t=0.16s
for i=1:1:t
  v(i)=w(0.4/dx+1,0.4/dy+1,i); %First, we have the temperature value
endfor
timevector=[0:dt:T-dt]
figure 2
plot(timevector,v);
xlabel('t (sec)')
ylabel('temp')
title('for x=y=0.4 and t=0.16s')
saveas(2,'x=y=0.4 and t=0.16s.png')

% Plotting the vertical temperature profile at t=0.16s and x=0.4
for j=1:1:Ny
  z(j)=w(j,0.4/dx+1,t); %First, we have the temperature value
endfor
figure 3
plot(x,z)
xlabel('x=0.4 (dimensionless)')
ylabel('temp')
title('Temperature at t=0.16s')
saveas(3,'x=0.4 and t=0.16s.png')
