clc
clear all
close all
%initializing dimensions for discretization in space
figure (1)
hold on
Nx = 2001;              %number of nodes in x
Ny = 81;                %number of nodes in y
len = 1;                %length
Re = 1e4;               %Reynold's number
h=0.1                   %height
dx = len/(Nx-1);
dy = h/(Ny-1);
L=1;
x_blasius=0.0005                   %x value chosen to compare with Blasius results

%initializing meshgrid

x = linspace(0,len,Nx);
y = linspace(0,h,Ny);
[X,Y]=meshgrid(x,y);

%% Initializing v and u vectors
v = zeros(Ny,Nx);
u = zeros(Ny,Nx);

%% Setting up boundary conditions
u(:,1) = 1; % u at inlet
v(:,1) = 0; % v at inlet
u(1,:) = 0; % No slip condition at the bottom
v(1,:) = 0; % No slip condition at the bottom
u(Ny,:) = 1; % Free outer flow at the top


%HERE I IMPLEMENT THE FDM

%% Numerical Solution using Finite Difference Method
for i = 1:Nx-1
    for j = 2:Ny-1  %determining u(i+1,j) from the derived FD equation
        u_sol(j) = u(j,i) + dx/(Re*u(j,i)*(dy^2))*(u(j+1,i)-2*u(j,i) +u(j-1,i)) - dx*v(j,i)/(2*dy*u(j,i))*(u(j+1,i)-u(j-1,i));
    end
    u(2:Ny-1,i+1) = u_sol(2:end); % determining v from the derived FD equation
    for j = 2:Ny
        v(j,i+1) = v(j-1,i+1) - dy/2/dx*(u(j,i+1)-u(j,i)+u(j-1,i+1)-u(j-1,i));
    end
end

%figure (1)
%contourf(x,y,u(:,:), 10)
%colorbar
%title('Contour plot for u')

figure (2)
contourf(x,y,v(:,:), 10)
colorbar
title('Contour plot for v')

%eta_little=y*100/(4*sqrt(0.005))
eta_big=y*sqrt(Re/L/x_blasius)
%figure (3)
%plot(y,u(:,2)*10)
%title('u velocity profile at x=0.0005')

figure (4)
plot(eta_big,u(:,(Nx-1)/2),'linewidth',2)
title('u velocity profile at x=0.0005 solving the boundary layer problem')
xlabel('\eta pour x = 0.0005', 'FontSize', 14);
ylabel(" u", 'FontSize', 14);
grid on;

%figure (5)
%plot(y,v(:,2)*10)
%title('v velocity profile at x=0.0005')

figure (6)
plot(eta_big,v(:,2),'linewidth',2)
title('v velocity profile at x=0.0005 solving the boundary layer problem')
xlabel('\eta pour x = 0.0005', 'FontSize', 14);
ylabel(" v", 'FontSize', 14);
grid on;
