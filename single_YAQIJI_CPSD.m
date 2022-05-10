clc;
clear;
close all;
addpath(genpath('E:\模态识别\基于非同步测量的压气机管道模态识别数据\Testcode-ver1\'));
addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\simulation\'));
addpath(genpath('E:\模态识别\模态识别代码\mode_detect\mode_detect\simulation\subprogram\'));
chemin = 'E:\模态识别\基于非同步测量的压气机管道模态识别数据\实验20-2019-11-16旋转机闸 测试';
%% 初始化参数%%%%%%%%%%%%% 
 
c=340; rho=1.225;pref=2e-5;
        zH = 0.4;         % 阵列的Z坐标
        NumMic = 12;       % 传声器的数量
        NumSM= 30;         % 非同步测量的次数
a=0.185; % 半径
S=pi*a^2; % 面积
Fs = 102400 ;      % 采样频率
time=5;            % 采样时间


%% %%%%%对测量信号进行处理%%%%%%%%%%%%

  nk=12;
  L_signal = Fs*time;
    L_seg = round(L_signal/100);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1); 
    rotor_speed=8000;
    
    eval(['load ''',chemin,'\','RotaryTest-10000-Rotate-No-1.mat''']); 
    Tdata_1=Data(:,1:13);
    
        Ind = [1:NumSM];   
         Num_file = Ind ;         
       for i_file =Num_file
    eval(['load ''',chemin,'\','RotaryTest-10000-Rotate-No-',num2str(i_file),'.mat''']);      
       Tdata=Data(:,1:13);
       phase_info{i_file}=cpsd(Tdata(:,13),Tdata_1(:,13),Wind,Noverlap,Nfft,Fs)./...
           abs(cpsd(Tdata(:,13),Tdata_1(:,13),Wind,Noverlap,Nfft,Fs));    
          for k=1:nk
            [temp,freq] = cpsd(Tdata(:,k),Tdata(:,13),Wind,Noverlap,Nfft,Fs);
              CC1(:,i_file+NumSM*(k-1)) = temp;
%            ./phase_info{i_file};
        end
 end
 %% 获取模态信息
 
    nk_enlarge=NumSM*nk;
    m=-nk_enlarge/2:nk_enlarge/2;
    for k=1:length(m)
        a_mf(:,k)=1/nk_enlarge*CC1(:,1:nk_enlarge)*exp(m(k)*1i*2*pi*(1:nk_enlarge)/nk_enlarge).'; 
    end
    %% 绘图

      h=figure('Visible', 'on');
%       set(gcf,'outerposition',get(0,'screensize'));%最大化
      set(gcf,'position',[200 100 800 600]);
     GAMMA =10*log10(abs(a_mf)/4e-10);
%      GAMMA(find(GAMMA<75))=50;
     imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA);ylim([1,rotor_speed/60*29*3.2]);
       xlim([-100,100]); 
     axis xy;