function [ay_floors,by_floors,cy_floors,model_prc_floors] = FloorsFacLoad_FRBAv_o(xt_tm1,vt_tm1,KAPPA,SIGMA,theta,rho0,rho1,rhov,lambda0,SIGMAlambda1,KAPPAv,SIGMAv,thetav,rho0_pi,rho1_pi,rhov_pi,sigq,gamv,gamvx,rho,Maturity,Strike)

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
%   ay = intercept vector (Nfloors x 1)
%   by = coefficient vector (Nmat x Nfac)
%   cy = coefficient vector (Nmat x 1)
%      where Nfloors=Nmat*Nk, ordered by "Maturity", then "Strike"
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


KAPPA_rn=KAPPA + SIGMAlambda1;
KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;
KAPPAv_rn=KAPPAv+sqrt(1-rho^2)*SIGMAv*gamv+rho*SIGMAv*gamvx;
[ay,by,cy] = YieldFacLoad_ODE3(Maturity,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA*SIGMA',rho0,rho1,rhov,KAPPAv_rn,SIGMAv^2);
ay=-ay./Maturity; ay=ay';
by=-by./repmat(Maturity,Nfac,1); by=by';
cy=-cy./Maturity; cy=cy';

ay_floors=zeros(Nfloors,1);
by_floors=zeros(Nfloors,Nfac);
cy_floors=zeros(Nfloors,1);
model_prc_floors=zeros(Nfloors,1);

for i_tau=1:Nmat
    tau=Maturity(i_tau);
    
    KAPPA_Q=KAPPA+SIGMAlambda1;
    theta_Q=inv(KAPPA_Q)*(KAPPA*theta-SIGMA*lambda0-SIGMA*SIGMA'*by(i_tau,:)'*tau);
    
    KAPPAv_Q=KAPPAv+sqrt(1-rho^2)*SIGMAv*gamv+rho*SIGMAv*gamvx+SIGMAv^2*cy(i_tau)*tau;
    thetav_Q=KAPPAv*thetav/KAPPAv_Q;
    
    rho0_Q=rho0_pi-sigq'*lambda0-sigq'*SIGMA'*by(i_tau,:)'*tau;
    rho1_Q=rho1_pi-lambda1'*sigq;
    rhov_Q=rhov_pi-gamvx-rho*SIGMAv*cy(i_tau)*tau;
        
    ax_Q=(eye(Nfac)-inv(KAPPA_Q*tau)*(eye(Nfac)-expm(-KAPPA_Q*tau)))*theta_Q;
    bx_Q=inv(KAPPA_Q*tau)*(eye(Nfac)-expm(-KAPPA_Q*tau));
    av_Q=(1-inv(KAPPAv_Q*tau)*(1-exp(-KAPPAv_Q*tau)))*thetav_Q;
    bv_Q=inv(KAPPAv_Q*tau)*(1-exp(-KAPPAv_Q*tau));  
    

    aI_Q=rho0_Q+rho1_Q'*ax_Q+rhov_Q*av_Q;
    bI_Q=rho1_Q'*bx_Q; bI_Q=bI_Q';
    cI_Q=rhov_Q*bv_Q;
    
    
    tmp_kron_Q=inv(kron(-KAPPA_Q,eye(Nfac))+kron(eye(Nfac),-KAPPA_Q));
    OMEGA_x_Q=reshape(tmp_kron_Q*reshape(expm(-KAPPA_Q*tau)*SIGMA*SIGMA'*expm(-KAPPA_Q'*tau)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

    sigq2=sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q;
    H0=sigq2'*sigq2-2*sigq2'*SIGMA'*inv(KAPPA_Q'*tau)*(eye(Nfac)-expm(-KAPPA_Q'*tau))*inv(KAPPA_Q')*rho1_Q  ...
        +rho1_Q'*inv(KAPPA_Q)*(1/tau*OMEGA_x_Q)*inv(KAPPA_Q')*rho1_Q;

    tmp1=rhov_Q*SIGMAv*(KAPPAv_Q)^(-1);
    H1=(1+2*rho*tmp1+tmp1^2)-2*(rho+tmp1)*tmp1*(KAPPAv_Q*tau)^(-1)*(1-exp(-KAPPAv_Q*tau)) ...
       +tmp1^2*(2*KAPPAv_Q*tau)^(-1)*(1-exp(-2*KAPPAv_Q*tau));

    H2=exp(-KAPPAv_Q*tau)*((1+2*rho*tmp1+tmp1^2)*(KAPPAv_Q*tau)^(-1)*(exp(KAPPAv_Q*tau)-1) ...
        -2*(rho+tmp1)*tmp1+tmp1^2*(KAPPAv_Q*tau)^(-1)*(1-exp(-KAPPAv_Q*tau)));
    
    dI_Q=H0+thetav_Q*(H1-H2);
    eI_Q=H2;


    for i_k=1:Nk
        strike_floors=Strike(i_k);

        tmp_mean=aI_Q+bI_Q'*xt_tm1+cI_Q*vt_tm1;
        tmp_var=dI_Q+eI_Q*vt_tm1;
        h0=tau*(tmp_mean+tmp_var/2);
        h1=(-log(1+strike_floors)+(tmp_mean+tmp_var))/sqrt(tmp_var/tau);
        h2=(-log(1+strike_floors)+tmp_mean)/sqrt(tmp_var/tau);

        h0x=tau*bI_Q; h1x=bI_Q/sqrt(tmp_var/tau); h2x=h1x;
        h0v=tau*(cI_Q+eI_Q/2);
        h1v=(cI_Q+eI_Q)/sqrt(tmp_var/tau)-1/2*eI_Q/tmp_var*h1;
        h2v=cI_Q/sqrt(tmp_var/tau)-1/2*eI_Q/tmp_var*h2;
        
        % floor pricing
        h1=-h1; h1x=-h1x; h1v=-h1v;
        h2=-h2; h2x=-h2x; h2v=-h2v;
        
        if ~(isreal(h1) & isreal(h2))
            warning('FLOORS pricing formula require real arguments');
            ay_floors=[];
            by_floors=[];
            cy_floors=[];
            return;
        end
        
        tmp_a=log(-exp(h0)*normcdf(h1)+(1+strike_floors)^tau*normcdf(h2));
        tmp_b=-(exp(h0)*(normcdf(h1)*h0x+normpdf(h1)*h1x)-(1+strike_floors)^tau*normpdf(h2)*h2x)/exp(tmp_a);        
        tmp_c=-(exp(h0)*(normcdf(h1)*h0v+normpdf(h1)*h1v)-(1+strike_floors)^tau*normpdf(h2)*h2v)/exp(tmp_a);

        id=i_k+(i_tau-1)*Nk;
        ay_floors(id)=-tau*ay(i_tau)+tmp_a-tmp_b'*xt_tm1-tmp_c*vt_tm1;
        by_floors(id,:)=-tau*by(i_tau,:)+tmp_b';
        cy_floors(id)=--tau*cy(i_tau)+tmp_c;
        
        model_prc_floors(id)=-tau*(ay(i_tau)+by(i_tau,:)*xt_tm1+cy(i_tau)*vt_tm1)+tmp_a;
    end
end








