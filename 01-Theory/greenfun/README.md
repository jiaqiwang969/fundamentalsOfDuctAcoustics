# Aim: reproduce the green function in the uniform flow duct
Ref: AIAA-2008-Rienstra-(20)、（17）



alpha_mn = Base.jmn_pm; 

kappa_p=(-w*M+sqrt(w^2-beta^2*alpha_mn.^2))/beta^2; 
kappa_m=(-w*M-sqrt(w^2-beta^2*alpha_mn.^2))/beta^2; 


Omega = w-kappa_mn;


Qp=+(kappa_p+Omega_p*M)*(1-m^2/alpha_p^2); 
Qm=-(kappa_p+Omega_m*M)*(1-m^2/alpha_m^2);


Gnp=-1/(2*pi*i)*bessel(m,alpha_p*r)*bessel(m,alpha_p*r0)/Qp/bessel(alpha_p)^2*exp(-i*kappa_p*(z-z0));
Gnm=-1/(2*pi*i)*bessel(m,alpha_m*r)*bessel(m,alpha_m*r0)/Qm/bessel(alpha_m)^2*exp(-i*kappa_m*(z-z0));


