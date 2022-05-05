% Fig5_deriv01.m
% Aim: prove Fig07
% Ref: Rienstra-(Fig07)
% 2021-05-05 wjq
%


clc
clear
close all

%% Add Subfunction
addpath(genpath('../chebfun-master'));
addpath(genpath('../subfunction'));


%% Sensor Positions
N = 100;
sensor.z0 = 8;
sensor.z = sensor.z0*ones(N,1);  % 
sensor.rho = ones(N,1); % all sensor in the wall
sensor.theta = 2*pi/N*[1:N].';
[sensor.x,sensor.y,sensor.z] = pol2cart(sensor.theta,sensor.rho,sensor.z);

%% Plot Duct Surface
rT = 1;h = 10;
[Xs,Ys,Zs] = cylinder(rT);
Zs = Zs*h;


%% Mode Generator
c = 343;               % sound speed
m = [0:22];              % circumferential modes maker
n = [1:10];
w = 20;                 % Non-dimensional frequency
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base] = BaseJ1(m,n(end),rT);  %Rienstra-50, at Least n-radialModes calculated
Eig = sqrt(w^2-Base.jmn_pm.^2) %Rienstra-52


figure
subplot(1,2,1)
hb=bar3(imag(Eig)); 
for j=1:length(hb)   
    zdata=get(hb(j),'Zdata');
    set(hb(j),'Cdata',zdata)
end
title("Imag(Eig)-right runing") 
colormap(jet);
colorbar  
% set(gca,'XTickLabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'},...
%     'yticklabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'})   
xlabel('m'); ylabel('n'); zlabel('z');  
box off  
grid on   
view([-90 90]);
subplot(1,2,2)
hb=bar3(real(Eig)); 
for j=1:length(hb)   
    zdata=get(hb(j),'Zdata');
    set(hb(j),'Cdata',zdata)
end
title("Re(Eig)-right runing") 
colormap(jet);
colorbar  
% set(gca,'XTickLabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'},...
%     'yticklabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'})   
xlabel('m'); ylabel('n'); zlabel('z');  
box off  
grid on   
view([-90 90]);


