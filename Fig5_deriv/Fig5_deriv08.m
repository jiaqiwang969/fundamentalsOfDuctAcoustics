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



% normalize
theta=0.87/3*pi %校准相位
a_mf_norm=a_mf./mean(abs(a_mf),2)*exp(i*theta);
scale=20;
a_mf_norm(find(abs(a_mf_norm)>scale))=scale;


m=[-nk_enlarge/2:nk_enlarge/2];
frequency=Sig.freq(1:15500);
amf=a_mf_norm(1:15500,:);

% 
% a_mf1=a_mf_norm;
% figure
% imagesc([-nk_enlarge/2:nk_enlarge/2],Sig.freq(1:15500),abs(a_mf1(1:15500,:)));...
%     ylim([1,Sig.rotor_speed/60*29*3.2]); axis xy;xlim([-100,100]);
% title(['截止模态分析',' -CPSD method '])
% colormap(jet);colorbar('peer',axes1);
% testTime='试验17-2019-11-11';
% title({[testTime,'-截止模态分析',' -旋转机匣测试 '];...
%     ['转速: ',num2str(Sig.rotor_speed),'rpm-采样率：',num2str(Sig.fs),'-',...
%     num2str(Sig.testPeriod),'s']},'FontSize',14)
% xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);
% set(axes1,'FontSize',14,'Layer','top');

figure
subplot(3,1,1)
imagesc(m,frequency,real(amf));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
colorbar
sgtitle(['exp-2019-11-11'],'FontSize',14)
title('real','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);

subplot(3,1,2)
imagesc(m,frequency,imag(amf));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
colorbar
title('imag','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);

subplot(3,1,3)
imagesc(m,frequency,abs(amf));
ylim([1,8000/60*29*3.2]); axis xy;xlim([-100,100]);
colormap(jet);
colorbar
title('imag','FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);



