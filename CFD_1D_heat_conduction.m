close all;
clear;
clc;

%Given Data
L = 1;                  %in m
alpha = 2.3*10^(-5);    %in m^2/sec
T_left = 100;             %in deg celcius
T_right = 200;          %in deg celcius

%Grid Parameters
nx = 101;               %number of grid points for x
nt = 20000;             %number of grid points for t
dx = L/(nx-1);          %in m
dt = 1;                 %in sec

%Initilizing the grid
T = zeros(nx,nt);

%Initial Condition
T(:,1) = 30;

%Boundary Conditions
T(1,:) = T_left;
T(end,:) = T_right;

error = 0;
epsilon = 1*10^(-4);    %error tolerance

%Using Crank-Nicholson Method

%Coefficients
a = -(alpha*dt)/(2*dt^2);
b = (1 + (alpha*dt)/(dt^2));
c = -(alpha*dt)/(2*dt^2);

%Tri diagonal Matrix
diagonals = [a*ones(nx-2, 1), b*ones(nx-2, 1), c*ones(nx-2, 1)];
offsets = [-1, 0, 1];    
A = spdiags(diagonals,offsets,nx-2,nx-2);
A
%K Matrix
K = zeros(nx-2,1);

for n = 1:nt-1
    T_new = T(:,n);
    for i = 2:nx-1
        K(i-1) = T(i,n) + alpha*dt*(T(i+1,n) - 2*T(i,n) + T(i-1,n))/dx^2;
    end

    T_new(2:nx-1) = A\K; 
    
    %error calculation
    error = max(abs(T_new - T(:,n+1)));
    
    %checking for convergence
    if error<epsilon
        break;
    end
    T(:,n+1) = T_new;
end


%Temperature profile with time evolution
figure(1)
x = linspace(0, L, nx);
for n = 1:1000:20000
    plot(x, T(:, n), 'DisplayName', sprintf('t = %d sec', n));hold on;
end

xlabel('Position (in m)');
ylabel('Temperature (in deg celcius)');
title('1-D Unsteady Heat Conduction');
legend('show');