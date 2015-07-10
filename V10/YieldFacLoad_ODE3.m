function [A,B,C] = YieldFacLoad_ODE3(TAUgrid,k,K,h,H,rho0,rho1,rhov,Kv,Hv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Solves the following system of ODEs
%   dA/dt = -rho0+k'*B+1/2*B'*H*B+h*C, A(0)=0,
%   dB/dt = -rho1-K'*B , B(0)=0,
%   dC/dt = -rhov-Kv*C+1/2*Hv*C^2, C(0)=0.
% -------------------------------------------------------
% INPUT:  
%   TAUgrid = ROW vector of maturities (1 x Ntau)
%   k       = COLUMN coefficient vector (N x 1)
%   K       = coefficient matrix (N x N)
% -------------------------------------------------------
% OUTPUT:
%   A = ROW solution vector (1 x Ntau)
%   B = MATRIX solution matrix (N x Ntau)
%   C = ROW solution matrix (1 x Ntau)
% -------------------------------------------------------
% NOTE: For the ith maturity in TAUgrid, A_const(i) and B_const(:,i) are
%   solution to the ODE system. 
% -------------------------------------------------------
% AUTHOR:   Bin Wei (Federal Resere Bank of Atlanta)
% DATE:     January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(k);
Ntau=length(TAUgrid);

[A1,B]=YieldFacLoad_ODE(TAUgrid,k,K,H,rho0,rho1);

% C satisfies: C'=(-rhov)+(-Kv)*C+(1/2*Hv)*C^2
%                =a+b*C+c*C^2
a=-rhov; b=-Kv; c=1/2*Hv;
if b^2-4*a*c>0
    ETA=sqrt(b^2-4*a*c);
else
    error('certain parameter restrictions are violated');
end
C=2*a*(1-exp(-ETA*TAUgrid))./((ETA+b)*exp(-ETA*TAUgrid)+(ETA-b));

if a==0
    A2=zeros(1,Ntau);
else
    A2=h*2*a*(TAUgrid/(ETA-b)+2/(ETA^2-b^2)*log(((ETA+b)*exp(-ETA*TAUgrid)+(ETA-b))/(2*ETA)));
end

A=A1+A2;