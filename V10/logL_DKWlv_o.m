function [logL,logL_vec,state,model_IE_options]= logL_DKWlv_o(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options)

BadLike= 9999;

pall=pallGet(pvar,paras_vec,pidx);
[T,Ny]=size(y_vec);
Nytips=length(TIPSgrid);
Noptions=length(MATgrid_options)*length(STRIKEgrid_options);

KAPPA=reshape(pall(1:Nfac^2),Nfac,Nfac);
SIGMA=reshape(pall(1+Nfac^2:2*Nfac^2),Nfac,Nfac);
theta=pall(2*Nfac^2+1:2*Nfac^2+Nfac);
rho0=pall(2*Nfac^2+Nfac+1);
rho1=pall(2*Nfac^2+Nfac+2:2*Nfac^2+2*Nfac+1);
rhov=pall(2*Nfac^2+2*Nfac+2);
lambda0=pall(2*Nfac^2+2*Nfac+3:2*Nfac^2+3*Nfac+2);
SIGMAlambda1=reshape(pall(2*Nfac^2+3*Nfac+3:3*Nfac^2+3*Nfac+2),Nfac,Nfac);
delta_y=pall(3*Nfac^2+3*Nfac+3:3*Nfac^2+3*Nfac+2+Ny);
delta_bcf=pall(3*Nfac^2+3*Nfac+3+Ny:3*Nfac^2+3*Nfac+4+Ny); % 2 forecast (6M and 1Y hor)
delta_bcfLT=pall(3*Nfac^2+3*Nfac+5+Ny);   

rho0_pi=pall(3*Nfac^2+3*Nfac+6+Ny);
rho1_pi=pall(3*Nfac^2+3*Nfac+7+Ny:3*Nfac^2+4*Nfac+6+Ny);
rhov_pi=pall(3*Nfac^2+4*Nfac+7+Ny);
sigq=pall(3*Nfac^2+4*Nfac+8+Ny:3*Nfac^2+5*Nfac+7+Ny);
KAPPAv=pall(3*Nfac^2+5*Nfac+8+Ny);
SIGMAv=pall(3*Nfac^2+5*Nfac+9+Ny);
thetav=pall(3*Nfac^2+5*Nfac+10+Ny);
gamv=pall(3*Nfac^2+5*Nfac+11+Ny);
gamvx=pall(3*Nfac^2+5*Nfac+12+Ny);
rho=pall(3*Nfac^2+5*Nfac+13+Ny);
delta_p=pall(3*Nfac^2+5*Nfac+14+Ny);
delta_tips=pall(3*Nfac^2+5*Nfac+15+Ny:3*Nfac^2+5*Nfac+14+Ny+Nytips);
delta_options=pall(3*Nfac^2+5*Nfac+15+Ny+Nytips:3*Nfac^2+5*Nfac+14+Ny+Nytips+length(MATgrid_options));

KAPPA_L=pall(3*Nfac^2+5*Nfac+15+Ny+Nytips+length(MATgrid_options));
SIGMA_L=pall(3*Nfac^2+5*Nfac+16+Ny+Nytips+length(MATgrid_options));
theta_L=pall(3*Nfac^2+5*Nfac+17+Ny+Nytips+length(MATgrid_options));
rho1_L=pall(3*Nfac^2+5*Nfac+18+Ny+Nytips+length(MATgrid_options):3*Nfac^2+6*Nfac+17+Ny+Nytips+length(MATgrid_options));
rhoL_L=pall(3*Nfac^2+6*Nfac+18+Ny+Nytips+length(MATgrid_options));
lambda0_L=pall(3*Nfac^2+6*Nfac+19+Ny+Nytips+length(MATgrid_options));
SIGMAlambda1_L=pall(3*Nfac^2+6*Nfac+20+Ny+Nytips+length(MATgrid_options));


lambda1=inv(SIGMA)*SIGMAlambda1;
lambda1_L=inv(SIGMA_L)*SIGMAlambda1_L;

if min(diag(KAPPA)) < 2.e-4 | KAPPAv <=0 | abs(rho)>1 | thetav<=1.e-12
    warning('KAPPA is too singular')
    logL=9999;
    logL_vec=[];
    state=[];
    model_prc_caps=[];
    return;
end

% Compute the nominal yield factor loadings
KAPPA_rn=KAPPA + SIGMAlambda1;
KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;
KAPPAv_rn=KAPPAv+sqrt(1-rho^2)*SIGMAv*gamv+rho*SIGMAv*gamvx;
if (-KAPPAv_rn)^2-4*(1/2*SIGMAv^2)*(-rhov)<=0  
     logL=987654321;
     logL_vec=[ ];
     state=[ ];
     model_IE_options=[];
     return;
end
[ay,by,cy] = YieldFacLoad_ODE3(MATgrid,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA*SIGMA',rho0,rho1,rhov,KAPPAv_rn,SIGMAv^2);
ay=-ay./MATgrid; ay=ay';
by=-by./repmat(MATgrid,Nfac,1); by=by';
cy=-cy./MATgrid; cy=cy';

% Compute real yield factor loadings
rho0_R=rho0-rho0_pi-1/2*sigq'*sigq+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
rhov_R=rhov-rhov_pi+gamvx-1/2;
lambda0_R=lambda0-sigq;
lambda1_R=lambda1;
SIGMAlambda1_R=SIGMA*lambda1_R;
KAPPA_rn_R=KAPPA + SIGMAlambda1_R;
KAPPAtheta_rn_R=KAPPA*theta - SIGMA*lambda0_R;
if (-KAPPAv_rn)^2-4*(1/2*SIGMAv^2)*(-rhov_R)<=0 
     logL=987654321;
     logL_vec=[ ];
     state=[ ];
     model_IE_options=[];
     return;
end
[ay_R,by_R,cy_R] = YieldFacLoad_ODE3(TIPSgrid,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA',rho0_R,rho1_R,rhov_R,KAPPAv_rn,SIGMAv^2);
ay_R=-ay_R./TIPSgrid; ay_R=ay_R';
by_R=-by_R./repmat(TIPSgrid,Nfac,1); by_R=by_R';
cy_R=-cy_R./TIPSgrid; cy_R=cy_R';


% Compute factor loadings for REAL component of TIPS yields
[Ay_TR,By_TR,Cy_TR]=YieldFacLoad_ODE3(TIPSgrid,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA',rho0_R,rho1_R+rho1_L,rhov_R,KAPPAv_rn,SIGMAv^2); 
ay_TR=-Ay_TR./TIPSgrid; ay_TR=ay_TR';
by_TR=-By_TR./repmat(TIPSgrid,Nfac,1); by_TR=by_TR';
cy_TR=-Cy_TR./TIPSgrid; cy_TR=cy_TR';


% Compute factor loadings for LIQUIDITY component of TIPS yields
KAPPA_L_rn=KAPPA_L + SIGMAlambda1_L;
KAPPAtheta_L_rn=KAPPA_L*theta_L - SIGMA_L*lambda0_L;
[Ay_L,By_L]=YieldFacLoad_ODE(TIPSgrid,KAPPAtheta_L_rn,KAPPA_L_rn,SIGMA_L*SIGMA_L',0,rhoL_L);
ay_L=-Ay_L./TIPSgrid; ay_L=ay_L';
by_L=-By_L./TIPSgrid; by_L=by_L';


ay3m=ay(1);   % first element should be 3M yield
by3m=by(1,:);
cy3m=cy(1);

ForecastHor=[0.5 1.0];
Nhor=length(ForecastHor);
af=zeros(Nhor,1);
bf=zeros(Nhor,Nfac);
cf=zeros(Nhor,1);

for i=1:Nhor
  af(i)=ay3m+by3m*(eye(Nfac)-expm(-KAPPA*ForecastHor(i)))*theta ...
         + cy3m*(1-exp(-KAPPAv*ForecastHor(i)))*thetav;
  bf(i,:)=by3m*expm(-KAPPA*ForecastHor(i));
  cf(i)=cy3m*exp(-KAPPAv*ForecastHor(i));
end

% z_{t+1}=A+B*z_t+eps^{z}_{t+1} where z_t=(p_t,x_t,v_t)
Bx=expm(-KAPPA*Dt);
Ax=(eye(Nfac)-Bx)*theta;

Bv=exp(-KAPPAv*Dt);
Av=(1-Bv)*thetav;

BL=exp(-KAPPA_L*Dt);
AL=(1-BL)*theta_L;

B=zeros(Nfac+3);
B(1,1)=1; B(1,2:(Nfac+1))=rho1_pi'*Dt; B(1,Nfac+2)=rhov_pi*Dt;
B(2:(Nfac+1),2:(Nfac+1))=Bx;
B(Nfac+2,Nfac+2)=Bv;
B(Nfac+3,Nfac+3)=BL;

A=zeros(Nfac+3,1);
A(1)=rho0_pi*Dt;
A(2:(Nfac+1))=Ax;
A(Nfac+2)=Av;
A(Nfac+3)=AL;


aI_Q=zeros(length(MATgrid_options),1);
bI_Q=zeros(length(MATgrid_options),Nfac);
cI_Q=zeros(length(MATgrid_options),1);

for i=1:length(MATgrid_options)
    MAT=MATgrid_options(i);

    [tmp_ay,tmp_by,tmp_cy] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA*SIGMA',rho0,rho1,rhov,KAPPAv_rn,SIGMAv^2);
    
    KAPPA_Q=KAPPA+SIGMAlambda1;
    theta_Q=inv(KAPPA_Q)*(KAPPA*theta-SIGMA*lambda0+SIGMA*SIGMA'*tmp_by);
    
    KAPPAv_Q=KAPPAv+sqrt(1-rho^2)*SIGMAv*gamv+rho*SIGMAv*gamvx-SIGMAv^2*tmp_cy;
    thetav_Q=KAPPAv*thetav/KAPPAv_Q;
    
    rho0_Q=rho0_pi-sigq'*lambda0+sigq'*SIGMA'*tmp_by;
    rho1_Q=rho1_pi-lambda1'*sigq;
    rhov_Q=rhov_pi-gamvx+rho*SIGMAv*tmp_cy;

        
    [tmp_aI_Q,tmp_bI_Q,tmp_cI_Q] = InfExpFacLoad_DKWv_option(KAPPA_Q,SIGMA,theta_Q,KAPPAv_Q,SIGMAv,thetav_Q,sigq,rho0_Q,rho1_Q,rhov_Q,rho,MAT);

    aI_Q(i)=tmp_aI_Q;
    bI_Q(i,:)=tmp_bI_Q; 
    cI_Q(i)=tmp_cI_Q; 
    
end


tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x_ss=-reshape(tmp_kron*reshape(SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);  
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Dt)*SIGMA*SIGMA'*expm(-KAPPA'*Dt)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

OMEGA_v_ss=thetav*SIGMAv^2/(2*KAPPAv);

OMEGA_L_ss=SIGMA_L^2/(2*KAPPA_L);
OMEGA_L=1/(2*KAPPA_L)*(1-exp(-2*KAPPA_L*Dt))*SIGMA_L^2;

OMEGA_q_ss=1e6;   % diffuse prior for price: choose a large enough number      
OMEGA_xq_ss=inv(KAPPA)*SIGMA*sigq;
OMEGA_vq_ss=rho*SIGMAv*thetav/KAPPAv;

OMEGA_ss=[OMEGA_q_ss, OMEGA_xq_ss',OMEGA_vq_ss,0;
       OMEGA_xq_ss,OMEGA_x_ss,zeros(Nfac,2);
       OMEGA_vq_ss,zeros(1,Nfac),OMEGA_v_ss,0;
       zeros(1,Nfac+2),OMEGA_L_ss];

IE_options_avg=[mean(IE_options(:,1:3)')',mean(IE_options(:,4:6)')'];

ik=find(p_vec(1:30));
priceinit=p_vec(ik(1));

% unconditional mean and variance
xtm1_tm1=[priceinit;theta;thetav;theta_L];  % E(z_t|t)
Ptm1_tm1=OMEGA_ss;                  % V(z_t|t)

logL_vec=zeros(T,1);
state=zeros(T,Nfac+3);
model_IE_options=zeros(T,length(MATgrid_options));

for j=1:T     
    
    % OMEGA: var[x_t|(t-1)]
    OMEGA_q=Dt*(sigq'*sigq + thetav)+(xtm1_tm1(Nfac+2)-thetav)*(1-exp(-KAPPAv*Dt))/KAPPAv;
    OMEGA_xq=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))*SIGMA*sigq;
    OMEGA_v=thetav*SIGMAv^2*(1-exp(-KAPPAv*Dt))^2/(2*KAPPAv)+xtm1_tm1(Nfac+2)*SIGMAv^2*(exp(-KAPPAv*Dt)-exp(-2*KAPPAv*Dt))/KAPPAv;
    OMEGA_vq=rho*SIGMAv*((1-exp(-KAPPAv/2*Dt))^2/KAPPAv*thetav+xtm1_tm1(Nfac+2)*(exp(-KAPPAv/2*Dt)-exp(-KAPPAv*Dt))/(KAPPAv/2));
    OMEGA=[OMEGA_q, OMEGA_xq',OMEGA_vq,0;
           OMEGA_xq,OMEGA_x,zeros(Nfac,2);
           OMEGA_vq,zeros(1,Nfac),OMEGA_v,0;
           zeros(1,Nfac+2),OMEGA_L];
       
    % prediction
    xt_tm1=A+B*xtm1_tm1;
    Pt_tm1=B*Ptm1_tm1*B'+OMEGA;

    % observation equation: y_t=a+b*z_t+eps^{y}_t
    if notmissing_p(j) == 0
        b= [zeros(Ny,1),by,cy,zeros(Ny,1)];
        a= ay;
        y_obs=y_vec(j,:)';
        Rvec=delta_y.^2;
    elseif notmissing_p(j) == 1
        b= [1,zeros(1,Nfac),0,0;zeros(Ny,1),by,cy,zeros(Ny,1)];
        a= [0;ay];
        y_obs=[p_vec(j); y_vec(j,:)'];
        Rvec=[delta_p^2; delta_y.^2];  
    end      

    if notmissing_bcf(j) == 1
        b = [b;[zeros(size(bf,1),1),bf,cf,zeros(size(bf,1),1)]];  
        a = [a;af];
        y_obs= [y_obs; bcf_vec(j,:)'];
        Rvec= [Rvec; delta_bcf.^2];
    end  
  
    if notmissing_bcfLT(j) == 1      
        hor=horLT(j);
        W=inv(KAPPA)*(expm(-KAPPA*hor)-expm(-KAPPA*(hor+5.0)))/5.0;
        Wv=inv(KAPPAv)*(exp(-KAPPAv*hor)-exp(-KAPPAv*(hor+5.0)))/5.0;
        afLT= ay3m + by3m*(eye(Nfac)-W)*theta+cy3m*(1-Wv)*thetav;
        bfLT= by3m*W;
        cfLT= cy3m*Wv;

        b = [b;[0,bfLT,cfLT,0]];
        a = [a;afLT];

        y_obs= [y_obs; bcfLT_vec(j)];
        Rvec= [Rvec; delta_bcfLT^2];
    end  

    if notmissing_tips(j) == 1
        b = [b;[zeros(size(by_TR,1),1),by_TR,cy_TR,by_L]];  
        a = [a;ay_TR+ay_L];
        y_obs= [y_obs; ytips_vec(j,:)'];
        Rvec= [Rvec; delta_tips.^2];
    end  
    
    if notmissing_options(j) == 1
        b = [b;[zeros(size(bI_Q,1),1),bI_Q,cI_Q,zeros(size(bI_Q,1),1)]];  
        a = [a;aI_Q];
        y_obs= [y_obs; IE_options_avg(j,:)'];
        Rvec= [Rvec; delta_options.^2];
        
        model_IE_options(j,:)=(aI_Q+bI_Q*xt_tm1(2:4)+cI_Q*xt_tm1(5))';
    end  
        

    R=diag(Rvec);

    yt_tm1=a+b*xt_tm1;
    Vt_tm1=b*Pt_tm1*b'+R;    
    
    et=y_obs-yt_tm1;    % forecast error
            
    rcVt=rcond(Vt_tm1);
    if (rcVt < 1.e-16 | isfinite(rcVt) == 0)
        warning('Vt is nearly singular')
        logL=9999;
        logL_vec=[ ];
        state=[ ];
        model_IE_options=[];
        return;
    end      

    logL_vec(j)=-0.5*(log(det(Vt_tm1))+et'*inv(Vt_tm1)*et);  
    
    % update the prediction of xt
    xtm1_tm1=xt_tm1+Pt_tm1*b'*inv(Vt_tm1)*et; 
    Ptm1_tm1=Pt_tm1-Pt_tm1*b'*inv(Vt_tm1)*b*Pt_tm1;
    
    state(j,:)=xtm1_tm1';
  
    if (xtm1_tm1(Nfac+2))<1e-12
        warning('volatility must be positive')
        logL=BadLike;
        logL_vec=[ ];
        state=[ ];
        model_IE_options=[];
        return;      
    end

end
logL= -(sum(logL_vec))/T;   % up to a constant term (add this)


if (isfinite(logL)==0 | imag(logL) ~=0 | isnan(logL))
  logL=9999;  
end

