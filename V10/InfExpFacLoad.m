function [aI,bI] = InfExpFacLoad(KAPPA,SIGMA,theta,sigq,sigqx,rho0_pi,rho1_pi,Maturity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Computes inflation expectation implied in inflation options
%   1/Maturity*log(E_t[Q_(t+Maturity)/Q_(t)])  = aI + bI'*x(t)  with
%   maturity contained in the "Maturity" vector. 
% -------------------------------------------------------
% INPUT:  
%   Maturity = ROW vector of maturities (1 x Nmat)
% -------------------------------------------------------
% OUTPUT:
%   aI = intercept vector (Nmat x 1)
%   bI = coefficient vector (Nfac x Nmat)
% -------------------------------------------------------
% AUTHOR:   Bin Wei (Federal Resere Bank of Atlanta)
% DATE:     May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nfac=length(theta);
tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Maturity)*SIGMA*SIGMA'*expm(-KAPPA'*Maturity)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);


bx=inv(KAPPA*Maturity)*(eye(Nfac)-expm(-KAPPA*Maturity));
ax=(eye(Nfac)-bx)*theta;

% 1/Maturity*Var_t{log(Q_(t+Maturity)/Q_(t))}
sigq2=sigq+SIGMA'*inv(KAPPA')*rho1_pi;
tmpIU=sigq2'*sigq2+sigqx^2 ...
    -2*sigq2'*SIGMA'*inv(KAPPA'*Maturity)*(eye(Nfac)-expm(-KAPPA'*Maturity))*inv(KAPPA')*rho1_pi  ...
    +rho1_pi'*inv(KAPPA)*(1/Maturity*OMEGA_x)*inv(KAPPA')*rho1_pi;

bI=bx'*rho1_pi;
aI=(rho0_pi+rho1_pi'*ax)+1/2*tmpIU; 

