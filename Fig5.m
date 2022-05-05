% Rienstar-DuctAcoustics-Fig5

clc
clear
addpath(genpath('chebfun-master'));
addpath(genpath('subfunction'));

frequency=[1000];
rT=1 %duct Radius
c=343;
m=[0]; %circumferential modes
Wave.k=5;%frequency*2*pi*rT/c-0.0000000001*i; %+-;  %Non-dimensional frequency

figure;
for n=[1 2 3 4 5 6]        %radial modes
    [Base]=BaseJ1(m,n,rT); %Rienstra-50
    Eig=sqrt(Wave.k^2-Base.jmn_pm.^2) %Rienstra-52
    plot(imag(i*Eig),real(i*Eig),'square');hold on
    plot(-imag(i*Eig),-real(i*Eig),'square');
end
xlim([-8 8])
ylim([-15 15])



% Duct Modes

%%
% mode=[-3:3];n_len=3;
%  N =131;Ratio=0.01; [D,r] = cheb(N,Ratio,1);
%  Mx=0.0*ones(N+1,1);%Mx=0.9-r.^2*0.9;%Figure11 (b)
% frequency=[1000];
% rT=1
% c=343;
% [initialEigValue,mode_enlarge,cutOffLine,len]=wm2initialEigValue(N,D,r,Ratio,Mx,Wave.k,Geo_b.m,Geo_b.n-1);
% frequency*2*pi*rT/c



