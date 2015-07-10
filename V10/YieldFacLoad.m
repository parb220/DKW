function [ay,by] = YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,Maturity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Computes intercept (A) and factor loadings (B) bonds with
%   maturity contained in the "Maturity" vector. 
% -------------------------------------------------------
% NOTE: A and B are such that P(t,n) = exp(A(n) + B(n)'*x(t))
%   where P(t,n) is the price at time t of a bond with maturity n, 
%   and x(t) is the state vector (Nfac x 1).
% -------------------------------------------------------
% NOTE: Affine Term Structure Model
%   - dx(t)=KAPPA*(theta-x(t))*dt+SIGMA*dB(t)
%   - r(t)=rho0+rho1'*x(t);
%   - the price of risk: lambda0+lambda1'*x(t)
%   where lambda1=inv(SIGMA)*SIGMAlambda1.
% -------------------------------------------------------
% INPUT:  
%   Maturity = ROW vector of maturities (1 x Nmat)
% -------------------------------------------------------
% OUTPUT:
%   ay = intercept vector (Nmat x 1)
%   by = coefficient vector (Nmat x Nfac)
% -------------------------------------------------------
% NOTE: Yield of zero-coupon bond with maturity n is given by
%   y(t,n)=ay(n)+by(n)*x(t)
% -------------------------------------------------------
% AUTHOR:   Bin Wei (Federal Resere Bank of Atlanta)
% DATE:     January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of state variables
Nfac=length(theta);

% Dynamics of the state vector under the risk-neutral measure
KAPPA_rn=KAPPA + SIGMAlambda1;
KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;

% A is a row vector (1 x Nmat)
% B is a matrix (Nfac x Nmat) 
[A,B]=YieldFacLoad_ODE(Maturity,KAPPAtheta_rn,KAPPA_rn,SIGMA*SIGMA',rho0,rho1);

% Convert from price to yield
ay=-A./Maturity;
by=-B./repmat(Maturity,Nfac,1);

ay=ay';
by=by';













