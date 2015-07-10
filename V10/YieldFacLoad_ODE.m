function [A,B] = YieldFacLoad_ODE(TAUgrid,k,K,H,rho0,rho1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Solves the following system of ODEs
%   dA/dt =  k'B + 1/2*B'H B - rho0 , A(0)=0,
%   dB/dt = -K'B - rho1 , B(0)=0.
% -------------------------------------------------------
% INPUT:  
%   TAUgrid = ROW vector of maturities (1 x Ntau)
%   k       = COLUMN coefficient vector (N x 1)
%   K       = coefficient matrix (N x N)
% -------------------------------------------------------
% OUTPUT:
%   A_const = ROW solution vector (1 x Ntau)
%   B_const = MATRIX solution matrix (N x Ntau)
% -------------------------------------------------------
% NOTE: For the ith maturity in TAUgrid, A_const(i) and B_const(:,i) are
%   solution to the ODE system. 
% -------------------------------------------------------
% AUTHOR:   Bin Wei (Federal Resere Bank of Atlanta)
% DATE:     January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(k);
Ntau=length(TAUgrid);

A = zeros(1,Ntau) ;
B = zeros(N,Ntau) ;

expk=cell(Ntau,1);
for i=1:Ntau
  expk{i,1}=expm(-K'*TAUgrid(i));
end

tmp=inv(-K')*(-rho1);
Kkron=inv(kron((-K),eye(N))+kron(eye(N),(-K)));
for i=1:Ntau
  B(:,i)=(expk{i,1}-eye(N))*tmp;
  tau=TAUgrid(i);
  K0=((expk{i,1}-eye(N))*inv(-K')-tau*eye(N))*tmp;
  K2=reshape(Kkron*reshape(expk{i,1}'*H*expk{i,1} - H,N^2,1),N,N);
  K1=H*inv(-K')*(expk{i,1}-eye(N));
  A(i)=k'*K0+0.5*tmp'*(K2-K1-K1'+tau*H)*tmp-tau*rho0;
  

%     % double check
%     [N_mat,D_mat]=eig(K);
%     Hstar=inv(N_mat)*H*inv(N_mat');
%     Hstar2=zeros(N,N);
%     for i=1:N
%         for j=1:N
%             Hstar2(i,j)=Hstar(i,j)*(1-exp(-(D_mat(i,i)+D_mat(j,j))*tau))/(D_mat(i,i)+D_mat(j,j));
%         end
%     end
%     K2=N_mat*Hstar2*N_mat';
  
end


