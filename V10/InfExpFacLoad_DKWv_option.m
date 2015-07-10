function [aI,bI,cI] = InfExpFacLoad_DKWv_option(KAPPA,SIGMA,theta,KAPPAv,SIGMAv,thetav,sigq,rho0_pi,rho1_pi,rhov_pi,rho,Maturity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Computes inflation expectation
%   1/Maturity*log(E_t[Q_(t+Maturity)/Q_(t)])=aI+bI'*x(t)+ cI*v(t) with
%   maturity contained in the "Maturity" vector. 
% -------------------------------------------------------
% NOTE: Affine Term Structure Model: q(t)=log(Q(t)) is log price level
%   - dx(t)=KAPPA*(theta-x(t))*dt+SIGMA*dB(t)
%   - dq(t)=pi(t)*dt+sigq'*dB(t)+sigqx*dW(t)
%   - pi(t)=rho0_pi+rho1_pi'*x(t)
% -------------------------------------------------------
% INPUT:  
%   Maturity = ROW vector of maturities (1 x Nmat)
% -------------------------------------------------------
% OUTPUT:
%   aI = intercept vector (Nmat x 1)
%   bI = coefficient vector (Nfac x Nmat )
% -------------------------------------------------------
% AUTHOR:   Bin Wei (Federal Resere Bank of Atlanta)
% DATE:     May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nfac=length(theta);
tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Maturity)*SIGMA*SIGMA'*expm(-KAPPA'*Maturity)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

bx=inv(KAPPA*Maturity)*(eye(Nfac)-expm(-KAPPA*Maturity));
ax=(eye(Nfac)-bx)*theta;

bv=inv(KAPPAv*Maturity)*(1-exp(-KAPPAv*Maturity));
av=(1-bv)*thetav;

sigq2=sigq+SIGMA'*inv(KAPPA')*rho1_pi;
H0=sigq2'*sigq2-2*sigq2'*SIGMA'*inv(KAPPA'*Maturity)*(eye(Nfac)-expm(-KAPPA'*Maturity))*inv(KAPPA')*rho1_pi  ...
    +rho1_pi'*inv(KAPPA)*(1/Maturity*OMEGA_x)*inv(KAPPA')*rho1_pi;

tmp1=rhov_pi*SIGMAv*(KAPPAv)^(-1);
H1=(1+2*rho*tmp1+tmp1^2)-2*(rho+tmp1)*tmp1*(KAPPAv*Maturity)^(-1)*(1-exp(-KAPPAv*Maturity)) ...
   +tmp1^2*(2*KAPPAv*Maturity)^(-1)*(1-exp(-2*KAPPAv*Maturity));

H2=exp(-KAPPAv*Maturity)*((1+2*rho*tmp1+tmp1^2)*(KAPPAv*Maturity)^(-1)*(exp(KAPPAv*Maturity)-1) ...
    -2*(rho+tmp1)*tmp1+tmp1^2*(KAPPAv*Maturity)^(-1)*(1-exp(-KAPPAv*Maturity)));

bI=bx'*rho1_pi;

cI=bv*rhov_pi+1/2*H2;

aI=(rho0_pi+rho1_pi'*ax+rhov_pi*av)+1/2*(H0+thetav*(H1-H2));



% tmp11=rho*rhov_pi*SIGMAv*(KAPPAv)^(-1);
% H11=(1+tmp11)^2*Maturity-2*(1+tmp11)*tmp11*(KAPPAv)^(-1)*(1-exp(-KAPPAv*Maturity)) ...
%     +tmp11^2*(2*KAPPAv)^(-1)*(1-exp(-2*KAPPAv*Maturity));
% 
% H12=exp(-KAPPAv*Maturity)*((1+tmp11)^2*(KAPPAv)^(-1)*(exp(KAPPAv*Maturity)-1) ...
%     -2*(1+tmp11)*tmp11+tmp11^2*(KAPPAv)^(-1)*(1-exp(-KAPPAv*Maturity)));
% 
% 
% tmp21=sqrt(1-rho^2)*rhov_pi*SIGMAv*(KAPPAv)^(-1);
% H21=tmp21^2*(Maturity-2*(KAPPAv)^(-1)*(1-exp(-KAPPAv*Maturity))+(2*KAPPAv)^(-1)*(1-exp(-2*KAPPAv*Maturity)));
% 
% H22=exp(-KAPPAv*Maturity)*tmp21^2*((KAPPAv)^(-1)*(exp(KAPPAv*Maturity)-1)-2*Maturity+(KAPPAv)^(-1)*(1-exp(-KAPPAv*Maturity)));
% 
% bI=1/Maturity*rho1_pi'*inv(KAPPA)*(eye(Nfac)-exp(-KAPPA*Maturity));
% 
% cI=1/Maturity*rhov_pi*(KAPPAv)^(-1)*(1-exp(-KAPPAv*Maturity))+1/(2*Maturity)*(H12+H22);
% 
% aI=(rho0_pi+rho1_pi'*theta+rhov_pi*thetav)-bI*theta-cI*thetav ...
%     +1/(2*Maturity)*(H0+thetav*(H11+H21));



% Alternatively way to calculate conditional variance: OMEGA_x
% [N_mat,D_mat]=eig(KAPPA);
% SIGMA2=inv(N_mat)*SIGMA;
% theta2=inv(N_mat)*theta;
% G0_mat=SIGMA2*SIGMA2';
% OMEGA2=zeros(Nfac,Nfac);
% OMEGA2_ss=zeros(Nfac,Nfac);
% 
% for i=1:Nfac
%     for j=1:Nfac
%         tmp=D_mat(i,i)+D_mat(j,j);
%         if (abs(tmp)<1e-20)
%             OMEGA2(i,j)=G0_mat(i,j)*Maturity;
%             OMEGA2_ss(i,j)=G0_mat(i,j)*Maturity;
%         else
%             OMEGA2(i,j)=G0_mat(i,j)*(1-exp(-tmp*Maturity))/tmp;
%             OMEGA2_ss(i,j)=G0_mat(i,j)/tmp;
%         end
%     end
% end
% OMEGA=N_mat*OMEGA2*N_mat';
% OMEGA_ss=N_mat*OMEGA2_ss*N_mat';














