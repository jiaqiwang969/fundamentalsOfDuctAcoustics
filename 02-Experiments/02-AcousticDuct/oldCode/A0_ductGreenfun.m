% Fig5_deriv02.m
% Aim: reconstruct simple form of green function
% ref: https://github.com/jiaqiwang969/SwirlDuctAcoustics/blob/main/SwirlDuctNoise.m
% 2021-05-09 wjq
%
% Experiment setup
% 48 microphones, 16/16/16
% Z=[0.939,0.839,0.539], theta=2pi/16*[1:16], rT=0.3
% two sources:r,theta,z=[0.2,0,0],[0.26,0,0]
% Aim: construct Green's function
% Nomalized step: 
% real   vs  normalized
% x=ax,t=at/c0,p=rho_0*c0^2*p,rho=rho_0*rho,v=c_0*v,
% 2pi*f=w*c0/a,k=k/a,Z=rho_0*c0*Z,U=c0*M

clc
clear
close all

%% Add Subfunction
addpath(genpath('../../01-Theory/chebfun-master'));
addpath(genpath('../../01-Theory/subfunction'));

rT=0.3
r_pole=[0.2 0.26]/rT;             %Pulse point position - upstream and downstream control
x_pole=linspace(-5,5,20);x_pole1=x_pole(find(x_pole<=0));x_pole2=x_pole(find(x_pole>0));
Entropy=0;                        %0:constant entropy condtion;1-2:logatithmic entropy condtion
Boundary=[1]; Type={'Hard Wall';'Lined Outer Wall';'Lined Inner Wall';'Lined Outer&Inner Wall'};
z_t=1-2*sqrt(-1);z_h=1-2*sqrt(-1);beta=[0.3];
N =31;Ratio=0.001; [D,r] = cheb(N,Ratio,1);

f=[6000]
w=5;%2*pi*rT/343;
Mx=0.5;
m=[-5:5];
n=[1:6];


[Base] = BaseJ1(m,n(end));   %Rienstra-50, at Least n-radialModes calculated
beta=sqrt(1-Mx^2);
Eig1=(-w*Mx-sqrt(w^2-beta^2*Base.jmn_pm.^2))/beta^2; %Rienstra-83
Eig2=(-w*Mx+sqrt(w^2-beta^2*Base.jmn_pm.^2))/beta^2; %Rienstra-83






for nk=1:length(m)

[G_nm1,Tgm11,Tgm12,Tgm13]=greenfun_dipoleNoise(r,Boundary,m(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,Eig1,z_t,z_h,r_pole,x_pole1,0,45,90);%green calculation£¨1*length(x_pole) cell£©
[G_nm2,Tgm21,Tgm22,Tgm23]=greenfun_dipoleNoise(r,Boundary,m(nk),Ratio,w,Tr,Omag,Mx,c02,rou0,Eig2,z_t,z_h,r_pole,x_pole2,0,45,90);%green calculation£¨1*length(x_pole) cell£©
[GNk1,TGm11,TGm12,TGm13]=cheb_cumKxCell(G_nm1,Tgm11,Tgm12,Tgm13,Ratio,length(x_pole1),length(mode1));%Sum all of wavenumbers*length(x_pole)
[GNk2,TGm21,TGm22,TGm23]=cheb_cumKxCell(G_nm2,Tgm21,Tgm22,Tgm23,Ratio,length(x_pole2),length(mode2));%Sum all of wavenumbers*length(x_pole)
[Gw1{1,nk},Tm11{1,nk},Tm12{1,nk},Tm13{1,nk},Gwn1{nk},TGmn11{nk},TGmn12{nk},TGmn13{nk}]=greenfun_span2volume(r,GNk1,TGm11,TGm12,TGm13,m,nk,x_pole1,thetaNumber);
[Gw2{1,nk},Tm21{1,nk},Tm22{1,nk},Tm23{1,nk},Gwn2{nk},TGmn21{nk},TGmn22{nk},TGmn23{nk}]=greenfun_span2volume(r,GNk2,TGm21,TGm22,TGm23,m,nk,x_pole2,thetaNumber);



end



figure
plot(r,V{mode1(2)}(4,:))
hold on
plot(r,U_r{1,1}(r))




figure
plot(lam(mode1),'r.');
hold on
plot(lam(mode2),'b.');
plot(Eig1,'r square')
plot(Eig2,'b square')

















% 
% 
% 
% 
% %% Sensor Positions
% N = 100;
% sensor.z0 = 8;
% sensor.z = sensor.z0*ones(N,1);  % 
% sensor.rho = ones(N,1); % all sensor in the wall
% sensor.theta = 2*pi/N*[1:N].';
% [sensor.x,sensor.y,sensor.z] = pol2cart(sensor.theta,sensor.rho,sensor.z);
% 
% %% Plot Duct Surface
% rT = 1;h = 10;
% [Xs,Ys,Zs] = cylinder(rT);
% Zs = Zs*h;
% h0=figure
% h0.Position = [100 400 540 400];
% surf(Zs,Xs,Ys,'MeshStyle','row',...
%     'MarkerFaceColor',[1 1 1],...
%     'FaceColor',[0.901960784313726 0.901960784313726 0.901960784313726],...
%     'EdgeColor',[0 0.447058823529412 0.741176470588235]);
% view([300 30]);
% grid('on');
% axis('tight');
% hold('on');
% axis equal
% xlabel({'x-axis'});
% ylabel({'y-axis'});
% xlabel({'z-axis'});
% title("Duct Acoustics")
% 
% %% Plot Sensor in the Duct
% scatter3(sensor.z,sensor.x,sensor.y,'o','filled','MarkerEdgeColor',[0 .5 .5],...
%     'MarkerFaceColor',[0 .7 .7],...
%     'LineWidth',15)
% 
% 
% %% Mode Generator
% c = 343;               % sound speed
% m = [-5:5];              % circumferential modes maker
% n = [1];
% w = 5;                 % Non-dimensional frequency
% % f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
% [Base] = BaseJ1(m,n(end),rT);  %Rienstra-50, at Least n-radialModes calculated
% Eig = sqrt(w^2-Base.jmn_pm.^2) %Rienstra-52
% 
% 
% 
% %% right runing modes with amplitude A_mn / left running modes B_mn
% % for simulation, it opens an "Interface" modified by yourself
% A_mn= rand(length(n),length(m))*2-1; A_mn= A_mn+i*(rand(length(n),length(m))*2-1);   
% B_mn= rand(length(n),length(m))*2-1; B_mn= B_mn+i*(rand(length(n),length(m))*2-1); 
% % A_mn=[
% % 0.117565484621756 + 0.467410026774460i
% % 0.147817414181606 + 0.656694003304474i
% % 0.0197646637288170 + 0.290185566632867i
% % 0.964291730299559 + 0.754536639716295i
% % 0.970372902013312 + 0.558118051358972i
% % 0.123860507146820 + 0.427792595548456i];
% % B_mn=[
% % 0.117565484621756 + 0.467410026774460i
% % 0.147817414181606 + 0.656694003304474i
% % 0.0197646637288170 + 0.290185566632867i
% % 0.964291730299559 + 0.754536639716295i
% % 0.970372902013312 + 0.558118051358972i
% % 0.123860507146820 + 0.427792595548456i];
% %A_mn(:,1)=0
% %B_mn(:,1)=0
% 
% %% plot amplitude
% h=figure
% h.Position = [100 100 540 200];
% h1=subplot(1,2,1)
% hb=bar3(abs(A_mn)); 
% for j=1:length(hb)   
%     zdata=get(hb(j),'Zdata');
%     set(hb(j),'Cdata',zdata)
% end
% title("right runing") 
% colormap("Parula");  
% colorbar  
% % set(gca,'XTickLabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'},...
% %     'yticklabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'})   
% xlabel('x'); ylabel('y'); zlabel('z');  
% box off  
% grid on   
% h1=subplot(1,2,2)
% hb=bar3(abs(B_mn)); 
% for j=1:length(hb)   
%     zdata=get(hb(j),'Zdata');
%     set(hb(j),'Cdata',zdata)
% end
% title("left runing") 
% colormap("Parula");  
% colorbar  
% % set(gca,'XTickLabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'},...
% %     'yticklabel',{'0.1',' ','0.3',' ','0.5',' ','0.7',' ','0.9'})   
% xlabel('x'); ylabel('y'); zlabel('z');  
% box off  
% grid on  
% 
% 
% % then adding those modes together
% % use "bsxfun" to speed up instead of by for-loop
% % for no-mean flow wave equation, Eig of right/left modes are \pm 
% 
% 
% % First, reproduce the solution of wave equation: %Rienstra-56(Normalized version)
% for km=1:length(m)
%     for kn=1:length(n)
%         U_r{kn,km}  = Base.normValue(n(kn),km)*chebfun(@(t) besselj(sign(m(km)).*m(km),Base.jmn_pm(n(kn))*t),[0,rT]);
%         pm{kn,km} =(A_mn(kn,km).*exp(-i*  Eig(n(kn),km)   *sensor.z0)+ ...
%                     B_mn(kn,km).*exp(-i* (-Eig(n(kn),km)) *sensor.z0)).*...
%                     U_r{kn,km}*exp(-i*m(km)*sensor.theta.'); % Rienstra-57a
%         um{kn,km} =(Eig(n(kn),km) * A_mn(kn,km).*exp(-i*   Eig(n(kn),km)   *sensor.z0)+...
%                    (Eig(n(kn),km))* B_mn(kn,km).*exp(-i* (-Eig(n(kn),km) )*sensor.z0)).*...
%                     U_r{kn,km}*exp(-i*m(km)*sensor.theta.')/w; % Rienstra-57b
% 
%     end
% end
% % p(z0,r,theta) & u(z0,r,theta) = \sum_{mn}{pm & um}
% pw = 0.0*pm{1,1};uw=0.0*um{1,1};
% for km = 1:length(m)
%     for kn=1:length(n)
%      pw = pw+pm{kn,km};
%      uw = uw+um{kn,km};
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Vertifying
% r = linspace(0,1,30).';
% [Theta,Rho] = meshgrid(sensor.theta,r');
% [yy,xx] = pol2cart(Theta,Rho);
% 
% 
% 
% h2=figure
% h2.Position = [650 400 440 400];
% offset = 0.05; cont = 22; % contour setting
% s1 = subplot(2,2,1); contour(xx,yy,real(pw(r)),cont); ...
%     axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
% title("pw") 
% s2 = subplot(2,2,2); contour(xx,yy,imag(pw(r)),cont); ...
%     axis('square'); xlabel(''); ylabel('imag', 'FontSize', 20);
% title("pw") 
% s1 = subplot(2,2,3); contour(xx,yy,real(uw(r)),cont); ...
%     axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
% title("uw") 
% s2 = subplot(2,2,4); contour(xx,yy,imag(uw(r)),cont); ...
%     axis('square'); xlabel(''); ylabel('imag', 'FontSize', 20);
% title("uw") 
% sgtitle(['m:[',num2str(m),'];n:[',num2str(n),']'], 'FontSize', 10);
% 
% 
% 
