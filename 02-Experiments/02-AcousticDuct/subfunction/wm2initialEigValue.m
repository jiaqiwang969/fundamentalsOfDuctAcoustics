function [initialEigValue,mode_enlarge,cutOffLine,len1]=wm2initialEigValue(N,D,r,Ratio,Mx,w,mode,n_len);
Entropy=0;Boundary=[1]; Type={'Hard Wall';'Lined Outer Wall';'Lined Inner Wall';'Lined Outer&Inner Wall'};
z_t=1-2*sqrt(-1);z_h=1-2*sqrt(-1);beta=[0.3]; 
Tr=0.0;Omag=0.0;
M_theta=Tr./r+Omag*r; %当旋流系数过大，则必须计算其影响;Mx不能为0
[c02,rou0,P0,s0]=entropyPara(r,N,Ratio,Omag,Tr,Entropy,beta);  
for wk=1:length(w)  
tic  
for nk=1:length(mode)
    [V{wk,nk},lam{wk,nk}]=eigfun_AB(r,D,N,w(wk),mode(nk),Ratio,Mx,M_theta,rou0,P0,c02,Boundary,z_t,z_h);%求特征
    crLayer=[min((w(wk)-mode(nk)*M_theta./r)./Mx);max((w(wk)-mode(nk)*M_theta./r)./Mx)];
    cutOffLine=GMM_Cluster3(lam{wk,nk},crLayer); 
    [mode1{wk,nk},len1{wk,nk}]=eig_choose_nmax(V{wk,nk},N,r,Ratio,lam{wk,nk},w(wk),Mx,Tr,Omag,-1,10100,-1,500,crLayer,cutOffLine,n_len,0); %选特征    
    disp(['m=',num2str(mode(nk)),'  w=',num2str(w(wk)),'  upstram/downstream   mode:  ',num2str(len1{wk,nk}),' kx:  ',num2str(lam{wk,nk}(mode1{wk,nk}).')]);   

end
toc
end

%% 调试代码段
fig = figure;
%subplot(1,2,1);
%title(['m=',num2str(mode)])
color1=hsv(length(mode));
for nk=1:length(mode)
for wk=1:length(w)
    %for mk=0:n_len
    %弄成三维的
    initialEigValue{wk,nk}=lam{wk,nk}(mode1{wk,nk});
    mode_enlarge{wk,nk}=mode(nk)*ones(size(lam{wk,nk}(mode1{wk,nk})));
    plot3(repmat(w(wk),length(mode1{wk,nk}),1),real(lam{wk,nk}(mode1{wk,nk})),repmat(mode(nk),length(mode1{wk,nk}),1),'.','Color',color1(nk,:));hold on;
    %plot3(repmat(w(wk),length(mode2{wk,nk}),1),real(lam{wk,nk}(mode2{wk,nk})),repmat(mode(nk),length(mode2{wk,nk}),1),'.','Color',color1(nk,:)); 
    %end
end
end 
hold on;
X = [0 0 w(end)  w(end);];
Y = [0 0 (-w(wk)*mean(Mx)./(1-mean(Mx).^2)) (-w(wk)*mean(Mx)./(1-mean(Mx).^2)); ];
Z = [mode(1) mode(end) mode(end)  mode(1); ];
para.TPP=fill3(X,Y,Z,'g','FaceAlpha',0.5); %添加颜色        % handle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set(para.TPP,'Visible','on')  
xlim([0 20]);ylim([-200 40]);grid on;xlabel('w'); ylabel('real(kx)');view(0,90); 
% dcm_obj = datacursormode(fig);datacursormode on,
% while (1)
%     presign = waitforbuttonpress;%如果按键盘为1，如果鼠标单击为0
%     k = 1;
%     if(presign==1)
%         c_info = getCursorInfo(dcm_obj);
%         a(k,:) = c_info.Position;
%         %通过c_info.Position识别特征值kx和特征函数V
%         i_w=find(w<c_info.Position(1)+0.001 & w>c_info.Position(1)-0.001);
%         i_m=find(mode==c_info.Position(3));
%         i_n=find(abs(real(lam{i_w,i_m}(mode1{i_w,i_m}))-c_info.Position(2))<0.0001);
%         disp(['m=',num2str( mode(i_m)),'  w=',num2str( w(i_w)),'  downstream mode:  ',num2str(len1{i_w,i_m}), ' kx:  ',num2str(lam{i_w,i_m}(mode1{i_w,i_m}).')]);   
%         subplot(1,2,2)
%         for ik=1:length(i_n)
%         plot(r,real(V{i_w,i_m}{mode1{i_w,i_m}(i_n(ik))}(1,:)));%V--》5*62，包含5个参数，分别是
%         hold on;
%         plot(r,imag(V{i_w,i_m}{mode1{i_w,i_m}(i_n(ik))}(1,:)),':');
%         end
%         iV1=real(V{i_w,i_m}{mode1{i_w,i_m}(i_n(ik))}(1,:))';iV1d=[0;diff(iV1)./diff(r)];iV1dd=diff(diff(diff(sign(iV1).*iV1d)));
%         iV2=imag(V{i_w,i_m}{mode1{i_w,i_m}(i_n(ik))}(1,:))';iV2d=[0;diff(iV2)./diff(r)];iV2dd=diff(diff(diff(sign(iV2).*iV2d)));   
%         iV=iV1-iV2;iVd=[0;diff(iV)./diff(r)];   iVdd=diff(diff(diff(sign(iV).*iVd)));
%         n_len=max(length(find(iVdd>abs(min(iVdd))&iVdd>0.5)),length(find(iV1dd>abs(min(iV1dd))&iV1dd>0.5)))
%          %n_len(kk)=max(length(find(iV1dd>abs(min(iV1dd))&iV1dd>0.5)),length(find(iV2dd>abs(min(iV2dd))& iV1dd>0.5)));%多次diff后，可以凸显root点
%          %n_len(kk)=max(length(find(diff(sign(iV1).*iV1d)<-1*max(diff(sign(iV1).*iV1d)))),length(find(diff(sign(iV2).*iV2d)<-1*max(diff(sign(iV2).*iV2d)))));%-4是根据经验选取（保证10阶径向模态正常）
%         tempMode=lam{i_w,i_m}(mode1{i_w,i_m});
%         
%         plot(r,iV1d);plot(r,[0;diff(sign(iV1).*iV1d)],'LineWidth',3);plot(r,[0;0;0;iVdd],'LineWidth',2);title(['轴向波数',num2str(tempMode(i_n)'),' 径向模态数：',num2str(n_len)]);
% %          figure;plot(r,imag(V{1,mode(kk)}(1,:)),'LineWidth',1.5);hold on;plot(r,iV2d);hold on;plot(r,sign(iV2).*iV2d);plot(r,[0;diff(sign(iV2).*iV2d)],'LineWidth',3);title(num2str(lam(mode(kk))));
% %          figure;plot(r,imag(V{1,mode(kk)}(1,:)),'LineWidth',1.5);hold on;plot(r,iVd);hold on;plot(r,sign(iV).*iVd);plot(r,[0;diff(sign(iV).*iVd)],'LineWidth',3);title(num2str(lam(mode(kk))));         
% %          figure;plot(iV1dd)
% %          figure;plot([0;diff(sign(iV2).*iV2d)])
%         hold off; 
%         k = k + 1;
%     end
% end

end









