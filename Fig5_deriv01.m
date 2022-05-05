% Fig5_deriv01.m
% Aim: restruct time signals for visual sensors
% 2021-05-05 wjq

clc
clear
close all

%% Add Subfunction
addpath(genpath('../chebfun-master'));
addpath(genpath('../subfunction'));


%% Sensor Positions
N = 7;
z = [8;8;8;8;8;8;8];
rho = [1;1;1;1;1;1;1];
theta = 2*pi/N*[0;1;2;3;4;5;6];
[x,y,z] = pol2cart(theta,rho,z);

%% Plot Duct Surface
rT = 1;h = 10;
[Xs,Ys,Zs] = cylinder(rT);
Zs = Zs*h;
surf(Zs,Xs,Ys,'MeshStyle','row',...
    'MarkerFaceColor',[1 1 1],...
    'FaceColor',[0.901960784313726 0.901960784313726 0.901960784313726],...
    'EdgeColor',[0 0.447058823529412 0.741176470588235]);
view([300 30]);
grid('on');
axis('tight');
hold('on');
axis equal
xlabel({'x-axis'});
ylabel({'y-axis'});
xlabel({'z-axis'});
title("Duct Acoustics")

%% Plot Sensor in the Duct
scatter3(z,x,y,'o','filled','MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',15)


%% Mode Generator

c = 343;               % sound speed
m = [4];               % circumferential modes maker
n = [3];
w = 5;                 % Non-dimensional frequency
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base]=BaseJ1(m,n,rT);       %Rienstra-50
Eig=sqrt(w^2-Base.jmn_pm.^2) %Rienstra-52


% Solution of wave equation: %Rienstra-54

for km=1:length(m)
    for kn=n
        phi_r=chebfun(@(t) besselj(sign(m(km)).*m(km),Base.jmn_pm(kn)*t),[0,rT]);
    end
end


%% Cut Plane of mode
zT =8;
t = 0:0.1:20;  % time
theta = linspace(0,2*pi,80);
r = linspace(0,1,30);

Phi_r = phi_r(r).';
psi=bsxfun(@times,exp(i*w*t),exp(-i*m*theta).')*exp(-i*Eig(n)*zT);
pm=bsxfun(@times,Phi_r,reshape(psi,1,length(theta),length(t)));
[Theta,Rho]=meshgrid(theta,r');
[yy,xx]=pol2cart(Theta,Rho);

figure
offset=0.05; cont=22; % contour setting
s1=subplot(1,2,1); contour(xx,yy,real(pm(:,:,4)),cont); ...
    axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
s2=subplot(1,2,2); contour(xx,yy,imag(pm(:,:,4)),cont); ...
    axis('square'); xlabel(''); ylabel('imag', 'FontSize', 20);
sgtitle(['Duct Mode-m',num2str(m),'-n',num2str(n)], 'FontSize', 30)
