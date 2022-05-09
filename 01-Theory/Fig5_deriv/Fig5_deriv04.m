% Fig5_deriv04.m
% Aim: restruct green function
% Ref: Rienstra-(120)
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
sensor.z0 = 0;
sensor.z = sensor.z0*ones(N,1);  % 
sensor.rho = ones(N,1); % all sensor in the wall
sensor.theta = 2*pi/N*[1:N].';
[sensor.x,sensor.y,sensor.z] = pol2cart(sensor.theta,sensor.rho,sensor.z);

%% Plot Duct Surface
rT = 1;h = 10;
[Xs,Ys,Zs] = cylinder(rT);
Zs = Zs*h;
h0=figure
h0.Position = [100 400 540 400];
surf(Zs,Xs,Ys,'MeshStyle','row',...
    'MarkerFaceColor',[1 1 1],...
    'FaceColor',[0.901960784313726 0.901960784313726 0.901960784313726],...
    'EdgeColor',[0 0.447058823529412 0.741176470588235]);
view([300 30]);
grid('on');
axis('tight');
hold('on');
axis equal
xlabel({'x-axis'});
ylabel({'y-axis'});
xlabel({'z-axis'});
title("Duct Acoustics")

%% Plot Sensor in the Duct
scatter3(sensor.z,sensor.x,sensor.y,'o','filled','MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',15)


%% Mode Generator
c = 343;               % sound speed
m = [-20:20];              % circumferential modes maker
n = [1:5];
w = 5;                 % Non-dimensional frequency
% f=[1000];            % physical frequency %f*2*pi*rT/c-0.0000000001*i;
[Base] = BaseJ1(m,n(end),rT);  %Rienstra-50, at Least n-radialModes calculated
Eig = sqrt(w^2-Base.jmn_pm.^2) %Rienstra-52




% First, reproduce the solution of wave equation: %Rienstra-56(Normalized version)
for km=1:length(m)
    for kn=1:length(n)
        temp1{kn,km}  = chebfun(@(t) besselj(sign(m(km)).*m(km),Base.jmn_pm(n(kn))*t),[0,rT]);
        temp2{kn,km} = temp1{kn,km}/(1-m(km)^2/Base.jmn_pm(n(kn))^2)/temp1{kn,km}(rT)/Eig(n(kn),km);
        pm{kn,km} =exp(-i*(-Eig(n(kn),km) )*abs(sensor.z0)).*temp2{kn,km}*exp(-i*m(km)*sensor.theta.')*w/(2*pi); % Rienstra-57a
    end
end
% p(z0,r,theta) & u(z0,r,theta) = \sum_{mn}{pm & um}
pw = 0.0*pm{1,1};
for km = 1:length(m)
    for kn=1:length(n)
     pw = pw+pm{kn,km};
    end
end










%% Vertifying
r = linspace(0,1,30).';
[Theta,Rho] = meshgrid(sensor.theta,r');
[yy,xx] = pol2cart(Theta,Rho);



h2=figure
h2.Position = [650 400 440 400];
offset = 0.05; cont = 22; % contour setting
% s1 = subplot(2,2,1); 
contour(xx,yy,real(pw(r)),cont); ...
    axis('square'); xlabel(''); ylabel('real', 'FontSize', 20);
title("pw") 



