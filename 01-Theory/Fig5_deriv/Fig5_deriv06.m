% Fig5_deriv05.m
% Aim: reconstruct several modes
% Ref: Rienstra-(8)
% 2021-05-05 wjq
% derived from Fig5_deriv02.m



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
m = [0:17];              % circumferential modes maker
n = [1:7];
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base] = BaseJ1(m,n(end),1);  %Rienstra-50, at Least n-radialModes calculated

w=linspace(0,20,1000);


figure
for kn=1:length(n)
    for k=1:length(w)
        Eig(k,:) = real(sqrt(w(k)^2-Base.jmn_pm(n(kn),:).^2))/rT; %Rienstra-52
    end

    for k=1:length(m)
        scatter3(w,Eig(:,k),n(kn)*ones(size(w)),'.');
        hold on
    end
end

ylabel({'Re(k_{m1})'});
xlabel({'w'});
title({'Figure-08'});
annotation('textbox',...
    [0.292857142857142 0.552380952380954 0.212946428571429 0.0547619047619048],...
    'String',{'Cut-on right running'});



