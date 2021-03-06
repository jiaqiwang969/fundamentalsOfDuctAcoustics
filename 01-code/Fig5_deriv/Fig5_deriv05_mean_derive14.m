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
m = [0:20];              % circumferential modes maker
n = [1];
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base] = BaseJ1(m,n(end),rT);  %Rienstra-50, at Least n-radialModes calculated
M=0.5; %mean flow
beta=sqrt(1-M^2);
%Eig = sqrt(w^2-Base.jmn_pm.^2) %Rienstra-52


w=linspace(0,20,10000);
for k=1:length(w)
%     Eig(k,:) = sqrt(w(k)^2-Base.jmn_pm.^2); %Rienstra-52
    Eigp(k,:)=(-w(k)*M+sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2; %Rienstra-83
    Eigm(k,:)=(-w(k)*M-sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2; %Rienstra-83

end



figure
for k=1:length(m)
   scatter(w,real(Eigp(:,k)),'.'); 
   hold on 
    scatter(w,real(Eigm(:,k)),'.');
end
ylabel({'Re(k_{m1})'});
xlabel({'w'});
title({'Figure-14'});
ylim([-20,20])
view(axes1,[-12.3 13.0878832116789]);

% annotation('textbox',...
%     [0.292857142857142 0.552380952380954 0.212946428571429 0.0547619047619048],...
%     'String',{'Cut-on right running'});
% 
% 

