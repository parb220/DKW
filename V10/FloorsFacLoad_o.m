function [ay_floors,by_floors,model_prc_floors] = FloorsFacLoad_o(xt_tm1,KAPPA,SIGMA,theta,rho0,rho1, ...
lambda0,SIGMAlambda1,rho0_pi,rho1_pi,sigq,sigqx,Maturity,Strike)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Computes intercept (A) and factor loadings (B) real bonds
%   with maturity contained in the "Maturity" vector. 
% -------------------------------------------------------
% NOTE: A, B, C are such that P(t,n)=exp(A(n)+B(n)'*x(t)+C(n)*v(t))
%   where P(t,n) is the price at time t of a bond with maturity n, 
%   and x(t) is the state vector (Nfac x 1) and v(t) is volatility.
% -------------------------------------------------------
% INPUT:  
%   Maturity = ROW vector of maturities (1 x Nmat)
%   Strike   = ROW vector of strikes (1 x Nk)
% -------------------------------------------------------
% OUTPUT:
%   ay = intercept vector (Ncaps x 1)
%   by = coefficient vector (Nmat x Nfac)
%   cy = coefficient vector (Nmat x 1)
%      where Ncaps=Nmat*Nk, ordered by "Maturity", then "Strike"
% -------------------------------------------------------
% NOTE: Log Price of CAPS with maturity n and strike k is approximated by
%   log(CAPS(t,n,k)=ay_CAPS(n,k)+by_CAPS(n,k)*x(t)+cy_CAPS(n,k)*v(t)
% -------------------------------------------------------
% AUTHOR:   Bin Wei (Federal Resere Bank of Atlanta)
% DATE:     February 22, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nfac=length(theta);     % Number of state variables
Nmat=length(Maturity);  % Number of maturities
Nk=length(Strike);      % Number of strike prices
Nfloors=Nmat*Nk;

lambda1=inv(SIGMA)*SIGMAlambda1;

% parameters under risk-neutral measure
KAPPA_rn=KAPPA + SIGMAlambda1;
KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;
theta_rn=inv(KAPPA_rn)*KAPPAtheta_rn;

[ay,by]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,Maturity);

ay_floors=zeros(Nfloors,1);
by_floors=zeros(Nfloors,Nfac);

model_prc_floors=zeros(Nfloors,1);

for i_tau=1:Nmat
    tau=Maturity(i_tau);
    
    KAPPA_Q=KAPPA+SIGMAlambda1;
    theta_Q=inv(KAPPA_Q)*(KAPPA*theta-SIGMA*lambda0-SIGMA*SIGMA'*by(i_tau,:)'*tau);
           
    rho0_Q=rho0_pi-lambda0'*sigq-sigq'*SIGMA'*by(i_tau,:)'*tau;
    rho1_Q=rho1_pi-lambda1'*sigq;
 
   
    ax_Q=(eye(Nfac)-inv(KAPPA_Q*tau)*(eye(Nfac)-expm(-KAPPA_Q*tau)))*theta;
%     bx_Q=(eye(Nfac)-expm(-KAPPA_Q'*tau))*inv(KAPPA_Q'*tau);    
    bx_Q=inv(KAPPA_Q*tau)*(eye(Nfac)-expm(-KAPPA_Q*tau));


    tmp_kron_Q=inv(kron(-KAPPA_Q,eye(Nfac))+kron(eye(Nfac),-KAPPA_Q));
    OMEGA_x_Q=reshape(tmp_kron_Q*reshape(expm(-KAPPA_Q*tau)*SIGMA*SIGMA'*expm(-KAPPA_Q'*tau)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

    H0=((sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q)'*(sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q)+sigqx^2)*tau ...
        -2*(sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q)'*SIGMA'*inv(KAPPA_Q')*(eye(Nfac)-expm(-KAPPA_Q'*tau))*inv(KAPPA_Q')*rho1_Q  ...
        +rho1_Q'*inv(KAPPA_Q)*OMEGA_x_Q*inv(KAPPA_Q')*rho1_Q;

    aI_Q=rho0_Q+rho1_Q'*ax_Q;
    bI_Q=rho1_Q'*bx_Q; bI_Q=bI_Q';
    
    for i_k=1:Nk
        strike_floors=Strike(i_k);

        tmp_mean=aI_Q+bI_Q'*xt_tm1;
        tmp_var=H0/tau;
        h0=tau*(tmp_mean+tmp_var/2);
        h1=(-log(1+strike_floors)+(tmp_mean+tmp_var))/sqrt(tmp_var/tau);
        h2=(-log(1+strike_floors)+tmp_mean)/sqrt(tmp_var/tau);

        h0x=tau*bI_Q; h1x=bI_Q/sqrt(tmp_var/tau); h2x=h1x;
        
        
        % floor pricing
        h1=-h1; h1x=-h1x; 
        h2=-h2; h2x=-h2x; 
        
        if ~(isreal(h1) & isreal(h2))
            warning('FLOORS pricing formula require real arguments');
            ay_floors=[];
            by_floors=[];
            return;
        end
        
        tmp_a=log(-exp(h0)*normcdf(h1)+(1+strike_floors)^tau*normcdf(h2));
        tmp_b=-(exp(h0)*(normcdf(h1)*h0x+normpdf(h1)*h1x)-(1+strike_floors)^tau*normpdf(h2)*h2x)/exp(tmp_a);        

        id=i_k+(i_tau-1)*Nk;
        ay_floors(id)=-tau*ay(i_tau)+tmp_a-tmp_b'*xt_tm1;
        by_floors(id,:)=-tau*by(i_tau,:)+tmp_b';
        
        model_prc_floors(id)=-tau*(ay(i_tau)+by(i_tau,:)*xt_tm1)+tmp_a;
    end
end








