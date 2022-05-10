function [G_nm,Tgm1,Tgm2,Tgm3]=greenfun_dipoleNoise(r,m,Ratio,w,Mx,k_nm,r0,x0,Angle1,Angle2,Angle3)
%numerical green's function calculation ,2018-10-10,wjq
diAxis1=[0,-cos(Angle1/180*pi),-sin(Angle1/180*pi)];%Direction of dipole force source
diAxis2=[0,-cos(Angle2/180*pi),-sin(Angle2/180*pi)];%Direction of dipole force source
diAxis3=[0,-cos(Angle3/180*pi),-sin(Angle3/180*pi)];%Direction of dipole force source

%p_nm=[];
G_nm=[];Tgm1=[];Tgm2=[];Tgm3=[];
if length(k_nm)&length(x0)~=0 %Wave number may be empty for high circumferential modes
    for kk=1:length(k_nm)
        xr = chebfun('xr', [Ratio, 1]);
        Mx = chebfun([fliplr(Mx')]',[Ratio,1]);
        %% for pn&Gn
        Phi=1-k_nm(kk)/w*Mx;
        A_rk=(-w^2*Phi^2)*Phi^2*w^2; %(James-3.5.21)
        B_rk=w^2*Phi^2*((-w^2*Phi^2)*(1/xr)+diff(w^2*Phi^2));%(James-3.5.22)
        C_rkw2=(-w^2*Phi^2)^2*(Phi^2-(k_nm(kk)/w)^2-(m/w)^2/xr^2);%(James-3.1.21)
        L=chebop(Ratio, 1);
        L.op=@(xr,g) A_rk*diff(g,2)+B_rk*diff(g,1)-w^2*C_rkw2.*g;%(James-3.5.20)
        L1=L;L2=L;

        L1.rbc=@(g) [g-1,diff(g)];g1 =L1\0;
        L2.lbc=@(g) [g-1,diff(g)];g2 =L2\0;

        %figure;plot(real(g1),'.-');hold on;plot(real(g2),'.-');
        %
        %%
        % df/dk
        Phi_k=-Mx;
        A_rk_k=(-w^2*2*Phi*Phi_k)*Phi^2*w^2+w^2*2*Phi*Phi_k*(-w^2*Phi^2);
        B_rk_k=w^2*2*Phi*Phi_k*((-w^2*Phi^2)*(1/xr)+diff(w^2*Phi^2))...
            +w^2*Phi^2*(-w^2*2*Phi*Phi_k*(1/xr)+diff(w^2*2*Phi*Phi_k));
        C_rk_kw2=2*(-w^2*Phi^2)*(-w^2)*2*Phi*Phi_k*(Phi^2-(k_nm(kk)/w)^2-(m/w)^2/xr^2)...
            +(-w^2*Phi^2)^2*(2*Phi*Phi_k-2*(k_nm(kk)/w));

        b_k1=w^2*C_rk_kw2*g1-B_rk_k*(diff(g1))-A_rk_k*(diff(g1,2));
        b_k2=w^2*C_rk_kw2*g2-B_rk_k*(diff(g2))-A_rk_k*(diff(g2,2));
        Lk1=L;Lk2=L;
        Lk1.rbc=@(g) [g-0,diff(g)];gk1 =Lk1\b_k1;
        Lk2.lbc=@(g) [g-0,diff(g)];gk2 =Lk2\b_k2;
        % figure
        % plot(real(gk1),'.-');hold on;plot(real(gk2),'.-');
        J=(w^2*Phi^2)*w^2*Phi^2;
        W=g1*diff(g2)-g2*diff(g1);
        J_k=2*w^4*Phi^3*Phi_k+(w^2*Phi^2)*w^2*2*Phi*Phi_k;
        W_k=gk1*diff(g2)+g1*diff(gk2)-g2*diff(gk1)-gk2*diff(g1); %(James-3.5.28)
        WJ_k=W_k*J+W*J_k;
        Amk=k_nm(kk)*Mx-w;Dmk=Amk^2;%(Posson2012AAIA-19)
        Tmk_1=diAxis1(3)*Dmk^2*k_nm(kk)+diAxis1(2)*(Dmk^2/xr*m);%(Posson2012AAIA-30)
        Tmk_2=diAxis2(3)*Dmk^2*k_nm(kk)+diAxis2(2)*(Dmk^2/xr*m);%(Posson2012AAIA-30)
        Tmk_3=diAxis3(3)*Dmk^2*k_nm(kk)+diAxis3(2)*(Dmk^2/xr*m);%(Posson2012AAIA-30)

        kkk=1;
        for xk0=x0
            %p_nm{kk,kkk}=sign(xk0)*((xr<r0)*sqrt(-1)*w/2/pi/r0/W_k(r0)*exp(sqrt(-1)*k_nm(kk)*xk0)*g1(r0)*g2...
            % +(xr>r0)*2*sqrt(-1)*w/4/pi/r0/W_k(r0)*exp(sqrt(-1)*k_nm(kk)*xk0)*g1*g2(r0));%(James-3.5.9)
            G_nm{kk,kkk}=sign(xk0)*((xr<r0)*sqrt(-1)*w/2/pi/r0/WJ_k(r0)*exp(sqrt(-1)*k_nm(kk)*xk0)*g1(r0)*g2...
                +(xr>r0)*2*sqrt(-1)*w/4/pi/r0/WJ_k(r0)*exp(sqrt(-1)*k_nm(kk)*xk0)*g1*g2(r0));%(James-3.5.9)
            G_m{kk,kkk}=sign(xk0)*((xr<r0)*sqrt(-1)*w/2/pi/r0/WJ_k(r0)*g1(r0)*g2...
                +(xr>r0)*2*sqrt(-1)*w/4/pi/r0/WJ_k(r0)*g1*g2(r0));%(James-3.5.9)
            Tgm1{kk,kkk}=sqrt(-1)*(Tmk_1*G_m{kk,kkk})*exp(sqrt(-1)*k_nm(kk)*xk0); %(Posson2012AAIA-44b),关于变量r0的函数，xd0=0
            Tgm2{kk,kkk}=sqrt(-1)*(Tmk_2*G_m{kk,kkk})*exp(sqrt(-1)*k_nm(kk)*xk0); %(Posson2012AAIA-44b),关于变量r0的函数，xd0=0
            Tgm3{kk,kkk}=sqrt(-1)*(Tmk_3*G_m{kk,kkk})*exp(sqrt(-1)*k_nm(kk)*xk0); %(Posson2012AAIA-44b),关于变量r0的函数，xd0=0

            kkk=kkk+1;
        end


    end
end


