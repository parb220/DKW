function [logL,logL_vec,state]= logL_DKWl_s(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,yS_vec,notmissing_swap)        

BadLike= 9999;

pall=pallGet(pvar,paras_vec,pidx);
[T,Ny]=size(y_vec);
Nytips=length(TIPSgrid);

KAPPA=reshape(pall(1:Nfac^2),Nfac,Nfac);
SIGMA=reshape(pall(1+Nfac^2:2*Nfac^2),Nfac,Nfac);
theta=pall(2*Nfac^2+1:2*Nfac^2+Nfac);
rho0=pall(2*Nfac^2+Nfac+1);
rho1=pall(2*Nfac^2+Nfac+2:2*Nfac^2+2*Nfac+1);
lambda0=pall(2*Nfac^2+2*Nfac+2:2*Nfac^2+3*Nfac+1);
SIGMAlambda1=reshape(pall(2*Nfac^2+3*Nfac+2:3*Nfac^2+3*Nfac+1),Nfac,Nfac);
delta_y=pall(3*Nfac^2+3*Nfac+2:3*Nfac^2+3*Nfac+1+Ny);
delta_bcf=pall(3*Nfac^2+3*Nfac+2+Ny:3*Nfac^2+3*Nfac+3+Ny); % 2 forecast (6M and 1Y hor)
delta_bcfLT=pall(3*Nfac^2+3*Nfac+4+Ny);   

rho0_pi=pall(3*Nfac^2+3*Nfac+5+Ny);
rho1_pi=pall(3*Nfac^2+3*Nfac+6+Ny:3*Nfac^2+4*Nfac+5+Ny);
sigq=pall(3*Nfac^2+4*Nfac+6+Ny:3*Nfac^2+5*Nfac+5+Ny);
sigqx=pall(3*Nfac^2+5*Nfac+6+Ny);
delta_p=pall(3*Nfac^2+5*Nfac+7+Ny);
delta_tips=pall(3*Nfac^2+5*Nfac+8+Ny:3*Nfac^2+5*Nfac+7+Ny+Nytips);
delta_swap=pall(3*Nfac^2+5*Nfac+8+Ny+Nytips:3*Nfac^2+5*Nfac+7+Ny+2*Nytips);
KAPPA_L=pall(3*Nfac^2+5*Nfac+8+Ny+2*Nytips);
SIGMA_L=pall(3*Nfac^2+5*Nfac+9+Ny+2*Nytips);
theta_L=pall(3*Nfac^2+5*Nfac+10+Ny+2*Nytips);
rho1_L=pall(3*Nfac^2+5*Nfac+11+Ny+2*Nytips:3*Nfac^2+6*Nfac+10+Ny+2*Nytips);
rhoL_L=pall(3*Nfac^2+6*Nfac+11+Ny+2*Nytips);
lambda0_L=pall(3*Nfac^2+6*Nfac+12+Ny+2*Nytips);
SIGMAlambda1_L=pall(3*Nfac^2+6*Nfac+13+Ny+2*Nytips);

lambda1=inv(SIGMA)*SIGMAlambda1;
lambda1_L=inv(SIGMA_L)*SIGMAlambda1_L;

if min(diag(KAPPA)) < 2.e-4
    warning('KAPPA is too singular')
    logL=9999;
    logL_vec=[];
    state=[];
    return;
end

% Compute NOMINAL yield factor loadings
[ay,by]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MATgrid);

% Compute REAL yield factor loadings (for real yields implied by inflation swaps)
rho0_R=rho0-rho0_pi-1/2*(sigq'*sigq+sigqx^2)+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
lambda0_R=lambda0-sigq;
SIGMAlambda1_R=SIGMAlambda1;
KAPPA_R_rn=KAPPA + SIGMAlambda1_R;
KAPPAtheta_R_rn=KAPPA*theta - SIGMA*lambda0_R;
[Ay_R,By_R]=YieldFacLoad_ODE(TIPSgrid,KAPPAtheta_R_rn,KAPPA_R_rn,SIGMA*SIGMA',rho0_R,rho1_R);
ay_R=-Ay_R./TIPSgrid; ay_R=ay_R';
by_R=-By_R./repmat(TIPSgrid,Nfac,1); by_R=by_R';

% Compute factor loadings for REAL component of TIPS yields
[Ay_TR,By_TR]=YieldFacLoad_ODE(TIPSgrid,KAPPAtheta_R_rn,KAPPA_R_rn,SIGMA*SIGMA',rho0_R,rho1_R+rho1_L);
ay_TR=-Ay_TR./TIPSgrid; ay_TR=ay_TR';
by_TR=-By_TR./repmat(TIPSgrid,Nfac,1); by_TR=by_TR';

% Compute factor loadings for LIQUIDITY component of TIPS yields
KAPPA_L_rn=KAPPA_L + SIGMAlambda1_L;
KAPPAtheta_L_rn=KAPPA_L*theta_L - SIGMA_L*lambda0_L;
[Ay_L,By_L]=YieldFacLoad_ODE(TIPSgrid,KAPPAtheta_L_rn,KAPPA_L_rn,SIGMA_L*SIGMA_L',0,rhoL_L);
ay_L=-Ay_L./TIPSgrid; ay_L=ay_L';
by_L=-By_L./TIPSgrid; by_L=by_L';

ay3=ay(1);   % first element should be 3M yield
by3=by(1,:);

ForecastHor=[0.5 1.0];
Nhor=length(ForecastHor);
af=zeros(Nhor,1);
bf=zeros(Nhor,Nfac);

for i=1:Nhor
  af(i)=ay3+by3*(eye(Nfac)-expm(-KAPPA*ForecastHor(i)))*theta;
  bf(i,:)=by3*expm(-KAPPA*ForecastHor(i));
end


% x_{t+1}=Ax+Bx*x_t+eps^{x}_{t+1}
Bx=expm(-KAPPA*Dt);
Ax= (eye(Nfac)-Bx)*theta;

BL=exp(-KAPPA_L*Dt);
AL=(1-BL)*theta_L;


B=zeros(Nfac+2);
B(1,1)=1; B(1,2:(Nfac+1))=rho1_pi'*Dt; 
B(2:(Nfac+1),2:(Nfac+1))=Bx;
B(Nfac+2,Nfac+2)=BL;

A=zeros(Nfac+2,1);
A(1)=rho0_pi*Dt;
A(2:(Nfac+1))=Ax;
A(Nfac+2)=AL;


tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x_ss=-reshape(tmp_kron*reshape(SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);  
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Dt)*SIGMA*SIGMA'*expm(-KAPPA'*Dt)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

OMEGA_L_ss=SIGMA_L^2/(2*KAPPA_L);
OMEGA_L=1/(2*KAPPA_L)*(1-exp(-2*KAPPA_L*Dt))*SIGMA_L^2;


OMEGA_ss=zeros(Nfac+2);
OMEGA_ss(1,1)=1e6;  % diffuse prior for price: choose a large enough number
OMEGA_ss(2:(Nfac+1),2:(Nfac+1))=OMEGA_x_ss;
OMEGA_ss((Nfac+2),(Nfac+2))=OMEGA_L_ss;

% OMEGA: var[x_t|(t-1)]
OMEGA_q=Dt*(sigq'*sigq + sigqx^2);
OMEGA_xq=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))*SIGMA*sigq;
OMEGA=[OMEGA_q, OMEGA_xq',0;
       OMEGA_xq,OMEGA_x,zeros(Nfac,1);
       zeros(1,Nfac+1),OMEGA_L];

ik=find(p_vec(1:30));
priceinit=p_vec(ik(1));

% unconditional mean and variance
xtm1_tm1=[priceinit;theta;theta_L];         % E(xt|t)
Ptm1_tm1=OMEGA_ss;                  % V(xt|t)


[T,Ny]=size(y_vec);
logL_vec=zeros(T,1);
state=zeros(T,Nfac+2);

err=zeros(T,Ny+Nytips+Nhor+2);
for j=1:T  
    % prediction
    xt_tm1=A+B*xtm1_tm1;
    Pt_tm1=B*Ptm1_tm1*B'+OMEGA;

    % observation equation: y_t=a+b*x_t+eps^{y}_t
    if notmissing_p(j) == 0
        b= [zeros(Ny,1),by,zeros(Ny,1)];
        a= ay;
        y_obs=y_vec(j,:)';
        Rvec=delta_y.^2;
    elseif notmissing_p(j) == 1
        b= [1,zeros(1,Nfac+1);zeros(Ny,1),by,zeros(Ny,1)];
        a= [0;ay];
        y_obs=[p_vec(j); y_vec(j,:)'];
        Rvec=[delta_p^2; delta_y.^2];  
    end      
  
    if notmissing_bcf(j) == 1
        b = [b;[zeros(size(bf,1),1),bf,zeros(size(bf,1),1)]];  
        a = [a;af];
        y_obs= [y_obs; bcf_vec(j,:)'];
        Rvec= [Rvec; delta_bcf.^2];
    end  

    if notmissing_bcfLT(j) == 1      
        hor=horLT(j);
        W=inv(KAPPA)*(expm(-KAPPA*hor)-expm(-KAPPA*(hor+5.0)))/5.0;
        afLT= ay3 + by3*(eye(Nfac)-W)*theta;
        bfLT= by3*W;

        b = [b;[0,bfLT,0]];
        a = [a;afLT];

        y_obs= [y_obs; bcfLT_vec(j)];
        Rvec= [Rvec; delta_bcfLT^2];
    end  

    if notmissing_tips(j) == 1
        b = [b;[zeros(size(by_TR,1),1),by_TR,by_L]];  
        a = [a;ay_L+ay_TR];
        y_obs= [y_obs; ytips_vec(j,:)'];
        Rvec= [Rvec; delta_tips.^2];
    end  
    
    for kk=1:Nytips
        if notmissing_swap(j,kk)
            b = [b;[0,by_R(kk,:),0]];  
            a = [a;ay_R(kk)];
            y_obs= [y_obs; yS_vec(j,kk)];
            Rvec= [Rvec; delta_swap(kk)^2];            
        end
    end
    
    R=diag(Rvec);

    yt_tm1=a+b*xt_tm1;
    Vt_tm1=b*Pt_tm1*b'+R;    
    et=y_obs-yt_tm1;    % forecast error

    rcVt=rcond(Vt_tm1);
%     if (rcVt < 1.e-16 | isfinite(rcVt) == 0)
    if (rcVt < 1.e-14 | isfinite(rcVt) == 0)
        warning('Vt is nearly singular')
        logL=9999;
        logL_vec=[ ];
        state=[ ];
        return;
    end      

    logL_vec(j)=-0.5*(log(det(Vt_tm1))+et'*inv(Vt_tm1)*et);  
    
    % update the prediction of xt
    xtm1_tm1=xt_tm1+Pt_tm1*b'*inv(Vt_tm1)*et; 
    Ptm1_tm1=Pt_tm1-Pt_tm1*b'*inv(Vt_tm1)*b*Pt_tm1;
    
    state(j,:)=xtm1_tm1';
end
logL= -sum(logL_vec)/T;   % up to a constant term (add this)


if (isfinite(logL)==0 | imag(logL) ~=0 | isnan(logL))
    logL=9999;
end


