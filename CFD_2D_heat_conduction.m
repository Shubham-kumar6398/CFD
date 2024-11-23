close all;
clear;
clc;

%Given Data
Lx = 1;                 %in m
Ly = 1;                 %in m
alpha = 2.3*10^(-5);    %in m^2/sec
Tb = 100;               %in deg celcius

%Grid Parameters
nx = 100;               %number of grid points for x
ny = 100;               %number of grid points for y
nt = 100000;            %number of grid points for t
dx = Lx/(nx-1);         %in m
dy = Ly/(ny-1);         %in m
dt = 1;

%Initilizing the grid
T = zeros(nx,ny,nt);

%Initial Condition
T(:,:,1) = 0;

%Boundary Conditions
T(1,:,:) = Tb;
T(end,:,:) = Tb;
T(:,1,:) = Tb;
T(:,end,:) = Tb;

error = 0;
epsilon = 1*10^(-4);    %error tolerance

%Using Explicit Method
for k = 2:nt
    for i = 2:nx-1
        for j = 2:ny-1
            T(i,j,k) = T(i,j,k-1) + (alpha*dt)*(T(i+1,j,k-1) -2*T(i,j,k-1) + T(i-1,j,k-1))/(dx^2) + (alpha*dt)*(T(i,j+1,k-1) -2*T(i,j,k-1) + T(i,j-1,k-1))/(dy^2);
        end
    end

    %error calculation
    error = abs(sum(T(:,:,k)) - sum(T(:,:,k-1)));

    %checking for convergence
    if error<epsilon
        disp('Convergence at t:')
        disp(k);
        break;
    end
    fprintf('Time step: %d / %d\n', k, nt);
end

%Final Heatmap
figure(1)
heatmap(T(:,:,k), 'Colormap', hot, 'ColorbarVisible', 'on')
xlabel('x');
ylabel('y');
title('2-D Unsteady Heat Conduction');

%Creating gif for the temperature evolution
figure(2)
for i = 1:100:k
    Temp = T(:,:,i);
    h = heatmap(Temp, 'Colormap', hot, 'ColorbarVisible', 'on');
    xlabel('x');
    ylabel('y');
    title('2-D Unsteady Heat Conduction');
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);

    if i == 1
        imwrite(A,map,'heat_animation_explicit.gif','gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,'heat_animation_explicit.gif','gif','WriteMode','append','DelayTime',0.1);
    end
end