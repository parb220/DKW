function [ay,by] = NomFacLoad(Kappa,Sigma,theta,rho_0,rho, lambda_1,DKappa,Maturity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Computes the intercept (A) and factor loadings (B) for all bonds
%   with maturity contained in the "Maturity" vector.  A and B
%   are such that P(t.n) = exp(A(n) + B(n)'*x(t)), where
%   P(t,n) is the price at time t of a bond with maturity n, 
%   and x(t) is the state vector.
%
% INPUTS:
%  Kappa,Sigma,theta,rho_0,rho,lambda_1,DKappa: structural parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nfac=length(theta);
RNK = Kappa + DKappa;
RNk = Kappa*theta - Sigma*lambda_1;

H=Sigma*Sigma';

[A,B]=ODE_FacLoad(Maturity,RNk,RNK,H,rho_0,rho);

% yields are given by y_{t,n}=ay(n) + By(n)*x_t
ay= -A./Maturity;
by= -B./repmat(Maturity,Nfac,1);

ay=ay';
by=by';



function [A_const,B_const] = ODE_FacLoad(T_grid,k,K,H,rho_0,rho)
% input: "T_grid" is a row-vector of time points
% output: A_const and B_const correspond to the solution to the following ODE
% dA/dt =  k'B + 1/2*B'H B - rho_0 , A(0)=0,
% dB/dt = -K'B - rho   , B(0)=0.

N=length(k);
times=length(T_grid);

A_const = zeros(1,times) ;
B_const = zeros(N,times) ;

expk=cell(times,1);
for i=1:times
  expk{i,1}=expm(-K'*T_grid(i));
end

Ka=-K';
iKa=inv(Ka);
krho=iKa*(-rho);

KK=inv(kron(Ka',eye(N))+kron(eye(N),Ka'));
for i=1:times
  B_const(:,i)=(expk{i,1}-eye(N))*krho;
  ti=T_grid(i);
  K2i=reshape(KK*reshape(expk{i,1}'*H*expk{i,1} - H,N^2,1),N,N);
  K1i=H*iKa*(expk{i,1}-eye(N));
  A_const(i)=k'*(iKa*(expk{i,1}-eye(N))-ti*eye(N))*krho ...
            +0.5*krho'*(K2i-K1i-K1i'+ti*H)*krho-ti*rho_0;
end







