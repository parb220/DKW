function [logL,logL_vec,state,model_prc_oil]= logL_FRBAoil(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,oil_vec,MATgrid_oil,notmissing_oil)       

BadLike= 9999;

pall=pallGet(pvar,paras_vec,pidx);
[T,Ny]=size(y_vec);
Nytips=length(TIPSgrid);
Noil=size(MATgrid_oil,1);

KAPPA=reshape(pall(1:Nfac^2),Nfac,Nfac);
SIGMA=reshape(pall(1+Nfac^2:2*Nfac^2),Nfac,Nfac);
theta=pall(2*Nfac^2+1:2*Nfac^2+Nfac);
rho0=pall(2*Nfac^2+Nfac+1);
rho1=pall(2*Nfac^2+Nfac+2:2*Nfac^2+2*Nfac+1);
rhod=pall(2*Nfac^2+2*Nfac+2);
rhos=pall(2*Nfac^2+2*Nfac+3);
lambda0=pall(2*Nfac^2+2*Nfac+4:2*Nfac^2+3*Nfac+3);
SIGMAlambda1=reshape(pall(2*Nfac^2+3*Nfac+4:3*Nfac^2+3*Nfac+3),Nfac,Nfac);
delta_y=pall(3*Nfac^2+3*Nfac+4:3*Nfac^2+3*Nfac+3+Ny);
delta_bcf=pall(3*Nfac^2+3*Nfac+4+Ny:3*Nfac^2+3*Nfac+5+Ny); % 2 forecast (6M and 1Y hor)
delta_bcfLT=pall(3*Nfac^2+3*Nfac+6+Ny);   

rho0_pi=pall(3*Nfac^2+3*Nfac+7+Ny);
rho1_pi=pall(3*Nfac^2+3*Nfac+8+Ny:3*Nfac^2+4*Nfac+7+Ny);
rhod_pi=pall(3*Nfac^2+4*Nfac+8+Ny);
rhos_pi=pall(3*Nfac^2+4*Nfac+9+Ny);
sigq=pall(3*Nfac^2+4*Nfac+10+Ny:3*Nfac^2+5*Nfac+9+Ny);
sigqx=pall(3*Nfac^2+5*Nfac+10+Ny);

SIGMAs=pall(3*Nfac^2+5*Nfac+11+Ny);
KAPPAd=pall(3*Nfac^2+5*Nfac+12+Ny);
SIGMAd=pall(3*Nfac^2+5*Nfac+13+Ny);
thetad=pall(3*Nfac^2+5*Nfac+14+Ny);
lambda0_d=pall(3*Nfac^2+5*Nfac+15+Ny);
lambda1_d=pall(3*Nfac^2+5*Nfac+16+Ny);

delta_p=pall(3*Nfac^2+5*Nfac+17+Ny);
delta_tips=pall(3*Nfac^2+5*Nfac+18+Ny:3*Nfac^2+5*Nfac+17+Ny+Nytips);
delta_oil=pall(3*Nfac^2+5*Nfac+18+Ny+Nytips);

lambda1=inv(SIGMA)*SIGMAlambda1;

if min(diag(KAPPA)) < 2.e-4 
    warning('KAPPA is too singular')
    logL=9999;
    logL_vec=[];
    state=[];
    model_prc_oil=[];
    return;
end

% parameters for log (oil) spot price
phi0=rho0-1/2*SIGMAs^2+SIGMAs*lambda0_d;
phi1=rho1;
phid=rhod-1+SIGMAs*lambda1_d;
phis=rhos;

% parameters under risk-neutral measure
KAPPA_rn=KAPPA + SIGMAlambda1;
KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;

KAPPAd_rn=KAPPAd+SIGMAd*lambda1_d;
thetad_rn=(KAPPAd*thetad-SIGMAd*lambda0_d)/KAPPAd_rn;

phi0_rn=rho0-1/2*SIGMAs^2;
phi1_rn=rho1;
phid_rn=rhod-1;
phis_rn=rhos;


tmp_k=[KAPPAtheta_rn;KAPPAd_rn*thetad_rn;phi0_rn];
tmp_K=[-KAPPA_rn',zeros(Nfac,1),phi1_rn;zeros(1,Nfac),-KAPPAd_rn,phid_rn;zeros(1,Nfac),0,phis_rn];
tmp_K=(-tmp_K)';
tmp_H=zeros(Nfac+2,Nfac+2); 
tmp_H(1:Nfac,1:Nfac)=SIGMA*SIGMA'; tmp_H(Nfac+1,Nfac+1)=SIGMAd^2; tmp_H(Nfac+2,Nfac+2)=SIGMAs^2;
[ay_tran,by_tran] = YieldFacLoad_ODE(MATgrid,tmp_k,tmp_K,tmp_H,rho0,[rho1;rhod;rhos]);

ay=ay_tran; by=by_tran(1:Nfac,:); dy=by_tran(Nfac+1,:); ey=by_tran(Nfac+2,:);

ay=-ay./MATgrid; ay=ay';
by=-by./repmat(MATgrid,Nfac,1); by=by';
dy=-dy./MATgrid; dy=dy';
ey=-ey./MATgrid; ey=ey';

% Compute real yield factor loadings
rho0_R=rho0-rho0_pi-1/2*(sigq'*sigq+sigqx^2)+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
rhod_R=rhod-rhod_pi;
rhos_R=rhos-rhos_pi;

lambda0_R=lambda0-sigq;
lambda1_R=lambda1;
SIGMAlambda1_R=SIGMA*lambda1_R;
KAPPA_rn_R=KAPPA + SIGMAlambda1_R;
KAPPAtheta_rn_R=KAPPA*theta - SIGMA*lambda0_R;

tmp_k_R=[KAPPAtheta_rn_R;KAPPAd_rn*thetad_rn;phi0_rn];
tmp_K_R=[-KAPPA_rn_R',zeros(Nfac,1),phi1_rn;zeros(1,Nfac),-KAPPAd_rn,phid_rn;zeros(1,Nfac),0,phis_rn];
tmp_K_R=(-tmp_K_R)';
[ay_tran_R,by_tran_R] = YieldFacLoad_ODE(TIPSgrid,tmp_k_R,tmp_K_R,tmp_H,rho0_R,[rho1_R;rhod_R;rhos_R]);
ay_R=ay_tran_R; by_R=by_tran_R(1:Nfac,:); dy_R=by_tran_R(Nfac+1,:); ey_R=by_tran_R(Nfac+2,:);

ay_R=-ay_R./TIPSgrid; ay_R=ay_R';
by_R=-by_R./repmat(TIPSgrid,Nfac,1); by_R=by_R';
dy_R=-dy_R./TIPSgrid; dy_R=dy_R';
ey_R=-ey_R./TIPSgrid; ey_R=ey_R';

ay3m=ay(1);   % first element should be 3M yield
by3m=by(1,:);
dy3m=dy(1);
ey3m=ey(1);


ForecastHor=[0.5 1.0];
Nhor=length(ForecastHor);
af=zeros(Nhor,1);
bf=zeros(Nhor,Nfac);
df=zeros(Nhor,1);
ef=zeros(Nhor,1);

for i=1:Nhor
  af(i)=ay3m+by3m*(eye(Nfac)-expm(-KAPPA*ForecastHor(i)))*theta ...
         + dy3m*(1-exp(-KAPPAd*ForecastHor(i)))*thetad ...
         + ey3m*(1/phis*(exp(phis*ForecastHor(i))-1)*(phi0+phi1'*theta+phid*thetad)) ...
         - ey3m*phi1'*inv(KAPPA+phis*eye(Nfac))*(exp(phis*ForecastHor(i))*eye(Nfac)-expm(-KAPPA*ForecastHor(i)))*theta ...
         - ey3m*phid*inv(KAPPAd+phis)*(exp(phis*ForecastHor(i))-exp(-KAPPAd*ForecastHor(i)))*thetad;         
  bf(i,:)=by3m*expm(-KAPPA*ForecastHor(i)) ...
         +ey3m*phi1'*inv(KAPPA+phis*eye(Nfac))*(exp(phis*ForecastHor(i))*eye(Nfac)-expm(-KAPPA*ForecastHor(i)));
  df(i)=dy3m*exp(-KAPPAd*ForecastHor(i)) ...
         +ey3m*phid*inv(KAPPAd+phid)*(exp(phis*ForecastHor(i))-exp(-KAPPAd*ForecastHor(i)));
  ef(i)=ey3m*exp(phis*ForecastHor(i));
end


% z_{t+1}=A+B*z_t+eps^{z}_{t+1} where z_t=(p_t,x_t,v_t)
Bx=expm(-KAPPA*Dt);
Ax=(eye(Nfac)-Bx)*theta;

Bd=exp(-KAPPAd*Dt);
Ad=(1-Bd)*thetad;


B=zeros(Nfac+3);
B(1,1)=1; B(1,2:(Nfac+1))=rho1_pi'*Dt; B(1,Nfac+2)=rhod_pi*Dt; B(1,Nfac+3)=rhos_pi*Dt;
B(2:(Nfac+1),2:(Nfac+1))=Bx;
B(Nfac+2,Nfac+2)=Bd;
B(Nfac+3,2:(Nfac+1))=phi1'*Dt; B(Nfac+3,Nfac+2)=phid*Dt; B(Nfac+3,Nfac+3)=1+phis*Dt;

A=zeros(Nfac+3,1);
A(1)=rho0_pi*Dt;
A(2:(Nfac+1))=Ax;
A(Nfac+2)=Ad;
A(Nfac+3)=phi0*Dt;

tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x_ss=-reshape(tmp_kron*reshape(SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);  
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Dt)*SIGMA*SIGMA'*expm(-KAPPA'*Dt)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

OMEGA_d_ss=SIGMAd^2/(2*KAPPAd);
OMEGA_s_ss=-SIGMAs^2/(2*phis);
OMEGA_ds_ss=SIGMAd*SIGMAs/(KAPPAd-phis);

OMEGA_q_ss=1e6;   % diffuse prior for price: choose a large enough number      
OMEGA_xq_ss=inv(KAPPA)*SIGMA*sigq;


OMEGA_ss=[OMEGA_q_ss, OMEGA_xq_ss',0,0;
       OMEGA_xq_ss,OMEGA_x_ss,zeros(Nfac,2);
       zeros(1,Nfac+1),OMEGA_d_ss,OMEGA_ds_ss;
       zeros(1,Nfac+1),OMEGA_ds_ss,OMEGA_s_ss];

ik=find(p_vec(1:30));
priceinit=p_vec(ik(1));

% unconditional mean and variance
xtm1_tm1=[priceinit;theta;thetad;-(phi0+phi1'*theta+phid*thetad)/phis];  % E(z_t|t)
Ptm1_tm1=OMEGA_ss;                  % V(z_t|t)

logL_vec=zeros(T,1);
state=zeros(T,Nfac+3);
model_prc_oil=zeros(T,Noil);

for j=1:T     
    % OMEGA: var[x_t|(t-1)]
    OMEGA_q=(sigq'*sigq + sigqx^2)*Dt;
    OMEGA_xq=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))*SIGMA*sigq;
    OMEGA_d=1/(2*KAPPAd)*(1-exp(-2*KAPPAd*Dt))*SIGMAd^2;
    OMEGA_s=1/(2*phis)*(exp(2*phis*Dt)-1)*SIGMAs^2;
    OMEGA_ds=1/(KAPPAd-phis)*(1-exp(-(KAPPAd-phis)*Dt))*SIGMAs*SIGMAd;
    
    OMEGA=[OMEGA_q, OMEGA_xq',0,0;
           OMEGA_xq,OMEGA_x,zeros(Nfac,2);
           zeros(1,Nfac+1),OMEGA_d,OMEGA_ds;
           zeros(1,Nfac+1),OMEGA_ds,OMEGA_s];
       
       
    % prediction
    xt_tm1=A+B*xtm1_tm1;
    Pt_tm1=B*Ptm1_tm1*B'+OMEGA;

    % observation equation: y_t=a+b*z_t+eps^{y}_t
    if notmissing_p(j) == 0
        b= [zeros(Ny,1),by,dy,ey];
        a= ay;
        y_obs=y_vec(j,:)';
        Rvec=delta_y.^2;
    elseif notmissing_p(j) == 1
        b= [1,zeros(1,Nfac+1),0;zeros(Ny,1),by,dy,ey];
        a= [0;ay];
        y_obs=[p_vec(j); y_vec(j,:)'];
        Rvec=[delta_p^2; delta_y.^2];  
    end      

    if notmissing_bcf(j) == 1
        b = [b;[zeros(size(bf,1),1),bf,df,ef]];  
        a = [a;af];
        y_obs= [y_obs; bcf_vec(j,:)'];
        Rvec= [Rvec; delta_bcf.^2];
    end  
  
    if notmissing_bcfLT(j) == 1      
        hor=horLT(j);
        W=inv(KAPPA)*(expm(-KAPPA*hor)-expm(-KAPPA*(hor+5.0)))/5.0;
        Wd=inv(KAPPAd)*(exp(-KAPPAd*hor)-exp(-KAPPAd*(hor+5.0)))/5.0;
        Ws=inv(phis)*(exp(phis*(hor+5.0))-exp(phis*hor))/5.0;
        
        afLT= ay3m + by3m*(eye(Nfac)-W)*theta ...
             +ey3m*1/phis*(Ws-1)*(phi0+phi1'*theta+phid*thetad) ...
             -ey3m*phi1'*inv(KAPPA+phis*eye(Nfac))*(Ws*eye(Nfac)-W)*theta ...
             -ey3m*phid*inv(KAPPAd+phid)*(Ws-Wd)*thetad;
        bfLT= by3m*W+ey3m*phi1'*inv(KAPPA+phis*eye(Nfac))*(Ws*eye(Nfac)-W);
        dfLT= dy3m*Wd+ey3m*phid*inv(KAPPAd+phis)*(Ws-Wd);
        efLT= ey3m*Ws;
        

        b = [b;[0,bfLT,dfLT,efLT]];
        a = [a;afLT];

        y_obs= [y_obs; bcfLT_vec(j)];
        Rvec= [Rvec; delta_bcfLT^2];
    end  

    
    if notmissing_oil(j)==1          
        tmp_k=[KAPPAtheta_rn;KAPPAd_rn*thetad_rn;phi0_rn+SIGMAs^2];
        tmp_K=[-KAPPA_rn',zeros(Nfac,1),phi1_rn;zeros(1,Nfac),-KAPPAd_rn,phid_rn;zeros(1,Nfac),0,phis_rn];
        tmp_K=(-tmp_K)';
        tmp_H=zeros(Nfac+2,Nfac+2); 
        tmp_H(1:Nfac,1:Nfac)=SIGMA*SIGMA'; tmp_H(Nfac+1,Nfac+1)=SIGMAd^2; tmp_H(Nfac+2,Nfac+2)=SIGMAs^2;
        [ay_tran_oil,by_tran_oil] = YieldFacLoad_ODE(MATgrid_oil(:,j)',tmp_k,tmp_K,tmp_H,(phi0_rn+1/2*SIGMAs^2),[phi1_rn;phid_rn;phis_rn]);
        
        Ay_oil=ay_tran_oil; By_oil=by_tran_oil(1:Nfac,:); 
        Dy_oil=by_tran_oil(Nfac+1,:); Ey_oil=by_tran_oil(Nfac+2,:)+1;

        Ay_oil=Ay_oil'; By_oil=By_oil'; Dy_oil=Dy_oil'; Ey_oil=Ey_oil';
        tmp_prc_oil=Ay_oil+By_oil*xt_tm1(2:Nfac+1)+Dy_oil*xt_tm1(Nfac+2)+Ey_oil*xt_tm1(Nfac+3);
        model_prc_oil(j,:)=tmp_prc_oil';
        
        b = [b;[zeros(size(By_oil,1),1),By_oil,Dy_oil,Ey_oil]];  
        a = [a;Ay_oil];
        y_obs = [y_obs; oil_vec(j,:)'];
        Rvec = [Rvec; delta_oil^2*ones(Noil,1)];      
    end
    
    if notmissing_tips(j) == 1
        b = [b;[zeros(size(by_R,1),1),by_R,dy_R,ey_R]];  
        a = [a;ay_R];
        y_obs= [y_obs; ytips_vec(j,:)'];
        Rvec= [Rvec; delta_tips.^2];
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
        model_prc_oil=[];
        return;
    end      

    logL_vec(j)=-0.5*(log(det(Vt_tm1))+et'*inv(Vt_tm1)*et);  
    
    % update the prediction of xt
    xtm1_tm1=xt_tm1+Pt_tm1*b'*inv(Vt_tm1)*et; 
    Ptm1_tm1=Pt_tm1-Pt_tm1*b'*inv(Vt_tm1)*b*Pt_tm1;
    
    state(j,:)=xtm1_tm1';
  
end
logL= -(sum(logL_vec))/T;   % up to a constant term (add this)


if (isfinite(logL)==0 | imag(logL) ~=0 | isnan(logL))
  logL=9999;  
end

