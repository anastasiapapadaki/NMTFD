% The 2D crank_Nickolson
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
Nxy=Nx*Ny;
alpha=dt/dx^2;
Result=zeros(Nx,Ny,t);

% Space discretization
x = 0:dx:1;
y = 0:dy:1;

%Initial Condition
w=zeros(Nxy,t); % we start with the initial temperature of zero everywhere

% Boundary conditions

for l=1:Nx:Nxy % w([1 1+Nx 1+2Nx... 1+Nx(Ny-1)], t) = 1-y3
  w(l,:) = 1-x((l-1)/Nx+1)^3;
end

for r=Nx:Nx:Nxy % w([Nx 2Nx... Nx*Ny)], t) = 1-sin(p/2 *y)
  w(r,:) = 1-sin(pi/2 * x(r/Nx));
end

for i=1:1:Nx
    w(i,:) = 1; % w(1->Nx, t) = 1
end
for j=0:1:Nx-1 %w(Nx(Ny-1)->Nx*Ny, t) = 1
    w(end-j,:) = 0;
end

%Crank-Nickolson penta diagonal matrix

diagonals = [2*(1+2*alpha)*ones(Nx*Ny,1), -alpha*ones(Nx*Ny,4)];
A=spdiags(diagonals,[0 -1 1 -Nx Nx],Nx*Ny,Nx*Ny);
I=speye(Nx*Ny);

%Replacement with identity lines into the matrix for the unchanging boundary
%conditions

for i=1:1:Ny
    A(i,:) = I(i,:);
    A(Nxy-Nx+i,:) = I(Nxy-Nx+i,:);
    A(i*Nx,:) = I(i*Nx,:);
    A(1+(i-1)*Nx,:) = I(1+(i-1)*Nx,:);
end
% In order to maintain the boundary conditions every element of boundary is
% facing the equivalent line of the identity


% Discretization
for k=1:t
    b=ones(Nx,1);%boundary conditions (w(1->Nx, t) = 1)
    for i=1:1:Ny-2;

        %stocking the precedent results into the result matrix to avoid a
        %second double loop

        Result(i,:,k)=w((i-1)*Nx+1:i*Nx,k);
        Result(Ny-1,:,k)=w((Ny-2)*Nx+1:(Ny-1)*Nx,k);
        Result(Ny,:,k)=w((Ny-1)*Nx+1:Ny*Nx,k);
        %Even if we could directly implement this border in Result doing
        % it in the main loop is a good way to check for bugs
        leftside=[alpha*w(1+i*Nx:Nx*(i+1)-2,k) + (2-4*alpha)*w(2+i*Nx:Nx*(i+1)-1,k)+ alpha*w(3+i*Nx:Nx*(i+1),k)+alpha*w(2+(i-1)*Nx:Nx*i-1,k)+alpha*w(2+(i+1)*Nx:(i+2)*Nx-1,k)];
        b=[b;1-x(i)^3;leftside; sin(1-sin(pi/2 * x(i)))];
    end
    b=[b;zeros(Nx,1)];%boundary conditions w(Nx(Ny-1)->Nx*Ny, t) = 1
    w(1:Nx*Ny,k+1)=A\b;
end
%Stocking the value of the last iteration
for i=1:1:Ny
    Result(i,:,t)=w((i-1)*Nx+1:i*Nx,t);
end
% Plotting the temperature surface at the final time

figure (1)
contourf(x,y,Result(:,:,t), 10)
w(1:Nxy,2)
colorbar

%Plotting the time evolution of the temperature at x=y=0.4

 v(:)=w((0.4/dx+1)*Nx+0.4/dx,:);
 timevector=[0:dt:T]
 figure (2)
 plot(timevector,v)
 xlabel('t(sec)')
 ylabel('temperature(dimensionless)')
 title('for x=y=0.4')
 saveas(4,"x=y=0.4 dt= 0.01 .png")

%Plotting the vertical temperature profile at t=0.16s and x=0.4
z(:)=w(Nx*0.4/dx+1:(Nx)*0.4/dx+Nx,0.16/dt);
figure (3)
plot(x,z)
xlabel('x=0.4 (dimensionless)')
ylabel('temperature (dimensionless)')
title('Temperature at t=0.16s and for x=0.4')
saveas(5,'x=0.4 and t=0.16s dt=0.01.png')
