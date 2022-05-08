clc
clear
close all
%% Experimental Data
Sig = importdata('S4-8000rpm-102400hz-5s-14probes_output_mat_CC1.mat'); 
nk_enlarge=length(Sig.CC2)*12;
m=-nk_enlarge/2:nk_enlarge/2;
for k=1:length(m)
    a_mf(:,k)=1/nk_enlarge*Sig.CC1(:,1:nk_enlarge)*...
        exp(m(k)*i*2*pi*(1:nk_enlarge)/nk_enlarge).';
end


%% 优化调试策略：
% 策略1: cut-off区域的能量占比尽量小
% 函数求解方法: 区域能量累加

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
rT = 0.19;h = 10;
[Xs,Ys,Zs] = cylinder(rT);
Zs = Zs*h;
%% Mode Generator
c = 343;               % sound speed
% m = [-100:100];        % circumferential modes maker
n = [1];
M = 0.0;               % mean flow in the duct
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base] = BaseJ1(m,n(end),1);  %Rienstra-50, at Least n-radialModes calculated
frequency=Sig.freq(1:15500);
w=frequency*2*pi*rT/c;%-0.0000000001*i;

for k=1:length(w)
    Eig(k,:) = sqrt(w(k)^2-(Base.jmn_pm).^2)*rT; %Rienstra-52
%     Eigp(k,:) = (-w(k)*M + sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2; % Rienstra-83-right running
%     Eigm(k,:) = (-w(k)*M - sqrt(w(k)^2-beta^2*Base.jmn_pm.^2))/beta^2; % Rienstra-83-left running
end

%frequency=w*c/(2*pi*rT);

% figure
z0=10*rT;
amf_theory=(-exp(i.*Eig*z0));






%%
% normalize
% scale=20;
%za_mf_norm(find(abs(a_mf_norm)>scale))=scale;
% kk=1;
% for scale=linspace(0.5,1.5,30)
a_mf_norm=a_mf(1:15500,:)./mean(abs(a_mf(1:15500,:)),2);




%总能量
E=sum(sum(abs(a_mf_norm)));

%内围能量
E0_plot=a_mf_norm.*(abs(amf_theory));
E0=sum(sum(abs(E0_plot)));

%外围能量
E1_plot=a_mf_norm.*(1-abs(amf_theory));
E1=sum(sum(abs(E1_plot)));

%raw/内围(propagator)
P0_plot=abs(a_mf(1:15500,:))./abs(E0_plot);
P0=sum(sum(abs(P0_plot)));

%costfun
costfunction=E0/E1;
%越大越好，说明分离的越彻底
% kk=kk+1;
% end



%% 截断
figure
subplot(1,3,1)
imagesc(m,frequency,abs(E0_plot));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
title('内围截断','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);

subplot(1,3,2)
imagesc(m,frequency,abs(E1_plot));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
title('外围截断','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);

subplot(1,3,3)
imagesc(m,frequency,log(P0_plot));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
title('cutoff-effct in log','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);




