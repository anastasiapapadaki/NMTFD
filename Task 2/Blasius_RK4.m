clc; clear all; close all; 
format short

%% RK4 to solve a system of 1st order ODEs 
%df1/d(eta)=f2 
%df2/d(eta)=f3
%df3/d(eta)=-1/2(f1*f3)
%where f1,f2,f3 are f,f',f'' respectively

%% initial conditions at h1=0.2
h1=0.2;
y1(1)=0;
y2(1)=0;
y3(1)=0.3323;

%% Defining our function handles using RK4 at h1=0.2
f1h1 = @(eta, y1, y2, y3) y2;
f2h1 = @(eta, y1, y2, y3) y3;
f3h1 = @(eta, y1, y2, y3) -1/2*y1*y3;
eta = 0:h1:14;
for i = 1:(length(eta)-1)
a= h1.*[f1h1(eta(i), y1(i), y2(i), y3(i)), f2h1(eta(i), y1(i), y2(i), y3(i)), f3h1(eta(i), y1(i), y2(i), y3(i))];
b= h1.*[f1h1(eta(i), y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f2h1(eta(i)+h1/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f3h1(eta(i)+h1/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2)];
c= h1.*[f1h1(eta(i), y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f2h1(eta(i)+h1/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f3h1(eta(i)+h1/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2)];
d= h1.*[f1h1(eta(i), y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f2h1(eta(i)+h1, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f3h1(eta(i)+h1, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3))];
y1(i+1) = y1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
y2(i+1) = y2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
y3(i+1) = y3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
end

%% initial conditions at h2=0.05
h2=0.05;
g1(1)=0;
g2(1)=0;
g3(1)=0.3323;


%% Defining our function handles using RK4 at h2=0.05
f1h2 = @(eta2, g1, g2, g3) g2;
f2h2 = @(eta2, g1, g2, g3) g3;
f3h2 = @(eta2, g1, g2, g3) -1/2*g1*g3;
eta2 = 0:h2:15;
for i = 1:(length(eta2)-1)
a= h2.*[f1h2(eta2(i), g1(i), g2(i), g3(i)), f2h2(eta2(i), g1(i), g2(i), g3(i)), f3h2(eta2(i), g1(i), g2(i), g3(i))];
b= h2.*[f1h2(eta2(i), g1(i)+a(1)/2, g2(i)+a(2)/2, g3(i)+a(3)/2), f2h2(eta2(i)+h2/2, g1(i)+a(1)/2, g2(i)+a(2)/2, g3(i)+a(3)/2), f3h2(eta2(i)+h2/2, g1(i)+a(1)/2, g2(i)+a(2)/2, g3(i)+a(3)/2)];
c= h2.*[f1h2(eta2(i), g1(i)+b(1)/2, g2(i)+b(2)/2, g3(i)+b(3)/2), f2h2(eta2(i)+h2/2, g1(i)+b(1)/2, g2(i)+b(2)/2, g3(i)+b(3)/2), f3h2(eta2(i)+h2/2, g1(i)+b(1)/2, g2(i)+b(2)/2, g3(i)+b(3)/2)];
d= h2.*[f1h2(eta2(i), g1(i)+c(1), g2(i)+c(2), g3(i)+c(3)), f2h2(eta2(i)+h2, g1(i)+c(1), g2(i)+c(2), g3(i)+c(3)), f3h2(eta2(i)+h2, g1(i)+c(1), g2(i)+c(2), g3(i)+c(3))];
g1(i+1) = g1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
g2(i+1) = g2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
g3(i+1) = g3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
end

%% Plot f
%figure(1)
%plot(eta ,y1,'-.r*')
%hold on
%plot(eta2 ,g1, '-.b')
%ylim([0 3])
%title(" Graph of f ");
%xlabel('\eta', 'FontSize', 14);
%ylabel("f ", 'FontSize', 14);
%grid on
%legend ('\eta =0.2', '\eta = 0.05')
%hold off
% saveas(figure(1),'F:\Computational Engineering\1st sem\Numerical Methods in Thermo-Fluid Dynamics I\Exercises\theoritical\task 2.4\RK4_graph of f','png')

%% Plot f'
%figure(2)
%plot(eta,y2,'-.r*')
%hold on
%plot(eta2,g2,'-.b')
%ylim([0 3])
%title(" Graph of f'");
%xlabel('\eta', 'FontSize', 14);
%ylabel("f' ", 'FontSize', 14);
%grid on
%legend ('\eta =0.2', '\eta = 0.05')
%hold off
% saveas(figure(2),'F:\Computational Engineering\1st sem\Numerical Methods in Thermo-Fluid Dynamics I\Exercises\theoritical\task 2.4\RK4_graph of fp','png')

%% Plot f''
%figure(3)
%plot(eta,y3,'-.r*')
%hold on
%plot(eta2,g3,'-.b')
%ylim([0 3])
%title(" Graph of f''");
%xlabel('\eta', 'FontSize', 14);
%ylabel("f''", 'FontSize', 14);
%grid on
%legend ('\eta =0.2', '\eta = 0.05')
%hold off
% saveas(figure(3),'F:\Computational Engineering\1st sem\Numerical Methods in Thermo-Fluid Dynamics I\Exercises\theoritical\task 2.4\RK4_graph of fpp','png')

%% Plot horizontal velocity u* obtained from RK4 at h2=0.05
%u*=f'(eta)=g2
figure(4)
plot(eta2,g2,'LineWidth',2)
ylim([0 1.4])
title(" u Velocity profile at x=0.5 using Blasius eqaution ");
xlabel('\eta pour x=0.00005', 'FontSize', 14);
xticks([0:1:15]);
ylabel(" u*", 'FontSize', 14);
grid on
% saveas(figure(4),'F:\Computational Engineering\1st sem\Numerical Methods in Thermo-Fluid Dynamics I\Exercises\theoritical\task 2.4\RK4_graph of horizonntal velocity','png')

%% Plot vertical velocity v* obtained from RK4 at h2=0.05

u_inf=10
L=1
Re=1e4
x_blasius=0.0005
v=0.5*sqrt(L/Re/x_blasius)*((eta2.*g2)-g1);
figure(5)
plot(eta2,v,'LineWidth', 2)
title(" v velocity profile at x=0.0005 using blasius equation");
xlabel('\eta = 0.05', 'FontSize', 14);
ylabel(" v*", 'FontSize', 14);
xticks([0:1:15]);
grid on
% saveas(figure(5),'F:\Computational Engineering\1st sem\Numerical Methods in Thermo-Fluid Dynamics I\Exercises\theoritical\task 2.4\RK4_graph of vertical velocity','png')

%%
%x = 0:0.05:8
%plot(x, 0.0345.*x.^(0.5))

%%
%p1 = plot(eta2, g1); hold on;
%p2 = plot(eta2, g2); 
%p3 = plot(eta2, g3); 
%p4 = plot(eta2, v); hold off;
%h = [p1;p2;p3;p4];
%legend(h,'f','f1 or u','f2','v', 'Location', 'Northwest')