%%----------------------Pwer by SJTU. CopyLicense------------------------%%

clc;clear;close all
addpath(genpath(pwd));
warning('off')
% % //======保存图像至指定文件夹=============== //
save_directory = ['./','结果输出/','声场模态输出',date];  %存储文件夹
if ~exist(save_directory) mkdir(save_directory);else disp('文件夹存在！');end

%%
load('database/mics_loc.mat')
% r,theta,x
% 三排传声器

mics_x=mics_loc(:,1).*cos(mics_loc(:,2));mics_y=mics_loc(:,1).*sin(mics_loc(:,2));%Xnet
mics_z=mics_loc(:,3);
%%


figure;plot3(mics_x,mics_y,mics_z,'o')
axis equal
hold on
% 声源面
plot3(0.2,0,0,'*'); plot3([0.2 0],[0 0],[0 0],'-');
plot3(0.26,0,0,'*');plot3([0.26 0],[0 0],[0 0],'-');
plot3([0 0],[0 0],[0 1.5],'--');

grid on

%%
Data = importdata('database/test20210121_test1_1_6000_1_1.mat');
Data_1=Data(:,1:16);Data_2=Data(:,17:32);Data_3=Data(:,33:48);

Fs = 16384;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = length(Data);             % Length of signal
t = (0:L-1)*T;        % Time vector

%% propagator matrix
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));
frequency=[5000];
rT = 0.3                            %duct Radius
c=343;
m=[-7:7];                       %circumferential modes
Wave.k=frequency*2*pi*rT/c; %+-;  %Non-dimensional frequency
n=[1:6]


[Base] = BaseJ1(m,n(end));   %Rienstra-50, at Least n-radialModes calculated
beta=sqrt(1-Mx^2);
Eig1=(-w*Mx-sqrt(w^2-beta^2*Base.jmn_pm.^2))/beta^2; %k+
Eig2=(-w*Mx+sqrt(w^2-beta^2*Base.jmn_pm.^2))/beta^2; %k-
Nmn2=(1-m^2/(Base.jmn_pm*rT)^2)*abs(besselj(sign(m(km)).*m(km),Base.jmn_pm(n(kn))*rT))^2;



%保证ref在同一位置，对三排分别做周向模态分解
refZ=mics_loc(1,3);
ref=Data(:,1);

for k = 1:15
    for km = 1:15
        P1(k,km) = exp(-i*m(km)*mics_loc(k,2));
        tmp = Base.normValue(1,km)*chebfun(@(t) besselj(sign(m(km)).*m(km),Base.jmn_pm(km)*t/rT),[0,rT]);
        P1(k,km) = tmp(rT)* P1(k,km) ;
    end
end

for k = 1:15
    for km = 1:15
        P2(k,km) = exp(-i*m(km)*mics_loc(k+16,2));
        tmp = Base.normValue(1,km)*chebfun(@(t) besselj(sign(m(km)).*m(km),Base.jmn_pm(km)*t/rT),[0,rT]);
        P2(k,km) = tmp(rT)* P2(k,km) ;
    end
end

for k = 1:15
    for km = 1:15
        P3(k,km) = exp(-i*m(km)*mics_loc(k+16*2,2));
        tmp = Base.normValue(1,km)*chebfun(@(t) besselj(sign(m(km)).*m(km),Base.jmn_pm(km)*t/rT),[0,rT]);
        P3(k,km) = tmp(rT)* P3(k,km) ;
    end
end


%% Duct Mode Decomposition
% method1:fft
% data_fft = 2^floor(log2(length(ref)));
% data = Data(1:data_fft,:);
% the_freq = [0:data_fft/2.56 - 1]*Fs/data_fft;  %数据频域离散刻度
% data_freq = fft(data)*2/data_fft;
% data_freq = data_freq(1:data_fft/2.56,:);

% method2:cpsd
L_signal = length(ref);
L_seg = round(L_signal/10);
Wind = hamming(L_seg);
Noverlap = round(L_seg/2);
Nfft = 2^(ceil(log2(L_seg))+1);
for k=1:48
    [temp,freq] = cpsd(Data(:,k),ref,Wind,Noverlap,Nfft,Fs);
    CC1(:,k) = temp;
end

%% 
amf1 = (inv(P1)*CC1(:,1:15).').';
a_mf1_norm=amf1./mean(abs(amf1),2);
amf2 = (inv(P2)*CC1(:,17:31).').';
a_mf2_norm=amf2./mean(abs(amf2),2);
amf3 = (inv(P3)*CC1(:,33:47).').';
a_mf3_norm=amf3./mean(abs(amf3),2);

figure
subplot(1,3,1)
imagesc(m,freq,abs((a_mf1_norm)));
ylim([1,8000]); axis xy; xlim([-7,7]);
colormap(jet);
colorbar
title('第一排abs','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);
subplot(1,3,2)
imagesc(m,freq,abs((a_mf2_norm)));
ylim([1,8000]); axis xy; xlim([-7,7]);
colormap(jet);
colorbar
title('第二排abs','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);
subplot(1,3,3)
imagesc(m,freq,abs((a_mf3_norm)));
ylim([1,8000]); axis xy; xlim([-7,7]);
colormap(jet);
colorbar
title('第三排abs','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);



