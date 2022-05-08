%% Clustering Using Gaussian Mixture Models
function [cutOffLine]=GMM_Cluster3(lam,crLayer) %10 36 -2 35
mode= find(real(lam)>-50000&real(lam)<50000&imag(lam)>1&imag(lam)<5000&abs(lam).^2>1E-3);
mode(find(crLayer(1)-0.1<real(lam(mode))&real(lam(mode))<crLayer(2)+0.1&abs(imag(lam(mode))))<1)=[];
cutOffLine=mean(real(lam(mode)));

