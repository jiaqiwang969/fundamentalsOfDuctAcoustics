function [Mode,Len]=eig_choose_nmax(V,N,r,Ratio,lam,w,Mx,Tr,Omag,re0,re1,im0,im1,crLayer,cutOffLine,n_max,cut_off);  %10 36 -2 35
%ѡ����:����������������critical layer
n_len=[];

    cmode=find(real(lam)<crLayer(2)+0.1&real(lam)>crLayer(1)-0.1&imag(lam)<5&imag(lam)>-5);
    mode=[1:length(lam)];

    for ii=1:length(cmode)
    mode(mode==cmode(ii))=[];
    end
%     Y(find(abs(imag(lam(Mode)))>15))=[];%cut-off mode �ų�
%     Mode(find(abs(imag(lam(Mode)))>15))=[];
     mode(find(abs(imag(lam(mode)))>15))=[];%cut-off mode �ų�
        
%     figure
%     for kk=1:length(Mode)
%     plot(imag(V{1,Mode(kk)}(1,:)'));hold on
%     end

    for kk=1:length(mode)
%       n_len(kk)=length(roots(chebfun(imag(V{1,mode(kk)}(1,:)'))));%ֱ��root������ܶ�0��������ֵ������Ϊ��root�㣬��Ҫ��diffɸѡ
          iV1=real(V{1,mode(kk)}(1,:))';iV1d=[0;diff(iV1)./diff(r)];iV1dd=diff(diff(diff(sign(iV1).*iV1d)));
          %iV2=imag(V{1,mode(kk)}(1,:))';iV2d=[0;diff(iV2)./diff(r)];iV2dd=diff(diff(diff(sign(iV2).*iV2d)));   
          iV=real(V{1,mode(kk)}(1,:))'-imag(V{1,mode(kk)}(1,:))';iVd=[0;diff(iV)./diff(r)];   iVdd=diff(diff(diff(sign(iV).*iVd)));
          n_len(kk)=max(length(find(iVdd>0.2*abs(min(iVdd))&iVdd>0.5)),length(find(iV1dd>0.2*abs(min(iV1dd))&iV1dd>0.5)));%���diff�󣬿���͹��root��
         %n_len(kk)=max(length(find(iV1dd>abs(min(iV1dd))&iV1dd>0.5)),length(find(iV2dd>abs(min(iV2dd))& iV1dd>0.5)));%���diff�󣬿���͹��root��
         %n_len(kk)=max(length(find(diff(sign(iV1).*iV1d)<-1*max(diff(sign(iV1).*iV1d)))),length(find(diff(sign(iV2).*iV2d)<-1*max(diff(sign(iV2).*iV2d)))));%-4�Ǹ��ݾ���ѡȡ����֤10�׾���ģ̬������
%          lam(mode(kk))
%          figure;plot(r,imag(V{1,mode(kk)}(1,:)),'LineWidth',1.5);hold on;plot(r,iV1d);hold on;plot(r,sign(iV1).*iV1d);plot(r,[0;diff(sign(iV1).*iV1d)],'LineWidth',3);title(num2str(lam(mode(kk))));
%          figure;plot(r,imag(V{1,mode(kk)}(1,:)),'LineWidth',1.5);hold on;plot(r,iV2d);hold on;plot(r,sign(iV2).*iV2d);plot(r,[0;diff(sign(iV2).*iV2d)],'LineWidth',3);title(num2str(lam(mode(kk))));
%          figure;plot(r,imag(V{1,mode(kk)}(1,:)),'LineWidth',1.5);hold on;plot(r,iVd);hold on;plot(r,sign(iV).*iVd);plot(r,[0;diff(sign(iV).*iVd)],'LineWidth',3);title(num2str(lam(mode(kk))));         
%          figure;plot(iV1dd)
%          figure;plot([0;diff(sign(iV2).*iV2d)])
    end
    
    mode(find((n_len>n_max)==1))=[];%����ģ̬�Ľ����������
    n_len(find((n_len>n_max)==1))=[];%����ģ̬�Ľ����������
    %[Y,I]=sort(n_len);Mode=mode(I);%����
    
    %cut-offģ̬������ȡ    
    t1=round(real(lam(mode))/0.00001)*0.00001;
    [B, I] = unique(t1, 'first'); %B����ź�������飬��ɾȥ�ظ��I�����������ԭλ��
    t2=setdiff(1:numel(t1), I); %ɾȥ�����λ��
    tk=[];tk=find(abs(t1)>200);
    for k=1:length(t2)
    tk=[tk;find(t1==t1(t2(k)))]; 
    end
    mode_cutoff=mode(tk);n_len_cutoff=n_len(tk);
    if cut_off==1;
    [Y_cutoff,I_cutoff]=sort(real(lam(mode(tk)))+0.001*imag(lam(mode(tk))));Mode=mode_cutoff(I_cutoff);Len=n_len_cutoff(I_cutoff);%���� %+0.1*imag(lam(mode))���ǵ�cut-off��imag���ֵ�����Ӱ��
    elseif cut_off==0;
    [Y,I]=sort(real(lam(mode))+0.001*imag(lam(mode)));Mode=mode(I);Len=n_len(I);%���� %+0.1*imag(lam(mode))���ǵ�cut-off��imag���ֵ�����Ӱ��
    end
    %��Ҫ�޸�fix������m=-1  w=10.1  downstream mode:  2  1  2  3  3  2  1  0
 end

