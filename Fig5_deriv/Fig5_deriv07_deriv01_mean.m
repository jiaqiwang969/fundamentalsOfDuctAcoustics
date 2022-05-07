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
% real experimenatl data 

rT = 0.185;h = 10;
[Xs,Ys,Zs] = cylinder(rT);
Zs = Zs*h;

%% Mode Generator
c = 343;               % sound speed
m = [-100:100];              % circumferential modes maker
n = [1];
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base] = BaseJ1(m,n(end),rT);  %Rienstra-50, at Least n-radialModes calculated
M=0.2; %mean flow
beta=sqrt(1-M^2);

w=linspace(0,42.0222,10000);
for k=1:length(w)
        %Eig(k,:) = (sqrt(w(k)^2-Base.jmn_pm.^2)); %Rienstra-52
        Eigp(k,:)=(-w(k)*M+sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2; %Rienstra-83
        Eigm(k,:)=(-w(k)*M-sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2; %Rienstra-83

end


frequency=w*c/(2*pi*rT);



[mm,ww]=meshgrid(m,frequency);

% figure
z0=0.2;
amfp=(exp(i.*Eigp*z0));
amfm=(exp(i.*Eigm*z0));

% offset = 0.05; cont = 22; % contour setting
% % s1 = subplot(2,2,1); 
% contour(mm,ww,abs(exp(i.*Eig*z0)),cont); ...
%     axis('square'); xlabel('m'); ylabel('f/hz', 'FontSize', 10);
% title("mode decomposition-what we see with experimental data") 
% grid on



figure
subplot(3,1,1)
imagesc(m,frequency,real(amfp));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
sgtitle(['ductMode-z0=',num2str(z0),'-r=',num2str(rT),'-M=',num2str(M)],'FontSize',14)
title('real','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);

subplot(3,1,2)
imagesc(m,frequency,imag(amfp));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
title('imag','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);

subplot(3,1,3)
imagesc(m,frequency,abs(amfp));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
title('imag','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);




