function [logL,logL_vec,state,model_prc_caps,model_prc_floors]= logL_FRBAv_o(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    caps_vec,MATgrid_caps,STRIKEgrid_caps,notmissing_caps, ...
    floors_vec,MATgrid_floors,STRIKEgrid_floors,notmissing_floors)       

BadLike= 9999;

pall=pallGet(pvar,paras_vec,pidx);
[T,Ny]=size(y_vec);
Nytips=length(TIPSgrid);
Ncaps=length(MATgrid_caps)*length(STRIKEgrid_caps);
Nfloors=length(MATgrid_floors)*length(STRIKEgrid_floors);

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
delta_options=pall(3*Nfac^2+5*Nfac+15+Ny+Nytips);

lambda1=inv(SIGMA)*SIGMAlambda1;

if min(diag(KAPPA)) < 2.e-4  
    warning('KAPPA is too singular')
    logL=9999;
    logL_vec=[];
    state=[];
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
    model_prc_caps=[];
    model_prc_floors=[];
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
     model_prc_caps=[];
     model_prc_floors=[];
     return;
end
[ay_R,by_R,cy_R] = YieldFacLoad_ODE3(TIPSgrid,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA',rho0_R,rho1_R,rhov_R,KAPPAv_rn,SIGMAv^2);
ay_R=-ay_R./TIPSgrid; ay_R=ay_R';
by_R=-by_R./repmat(TIPSgrid,Nfac,1); by_R=by_R';
cy_R=-cy_R./TIPSgrid; cy_R=cy_R';

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

B=zeros(Nfac+2);
B(1,1)=1; B(1,2:(Nfac+1))=rho1_pi'*Dt; B(1,end)=rhov_pi*Dt;
B(2:(Nfac+1),2:(Nfac+1))=Bx;
B(Nfac+2,Nfac+2)=Bv;

A=zeros(Nfac+2,1);
A(1)=rho0_pi*Dt;
A(2:(Nfac+1))=Ax;
A(Nfac+2)=Av;

tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x_ss=-reshape(tmp_kron*reshape(SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);  
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Dt)*SIGMA*SIGMA'*expm(-KAPPA'*Dt)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

OMEGA_v_ss=thetav*SIGMAv^2/(2*KAPPAv);

OMEGA_q_ss=1e6;   % diffuse prior for price: choose a large enough number      
OMEGA_xq_ss=inv(KAPPA)*SIGMA*sigq;
OMEGA_vq_ss=rho*SIGMAv*thetav/KAPPAv;

OMEGA_ss=[OMEGA_q_ss, OMEGA_xq_ss',OMEGA_vq_ss;
       OMEGA_xq_ss,OMEGA_x_ss,zeros(Nfac,1);
       OMEGA_vq_ss,zeros(1,Nfac),OMEGA_v_ss];

ik=find(p_vec(1:30));
priceinit=p_vec(ik(1));

% unconditional mean and variance
xtm1_tm1=[priceinit;theta;thetav];  % E(z_t|t)
Ptm1_tm1=OMEGA_ss;                  % V(z_t|t)

logL_vec=zeros(T,1);
state=zeros(T,Nfac+2);
model_prc_caps=zeros(T,Ncaps);
model_prc_floors=zeros(T,Ncaps);

for j=1:T  
    OMEGA_q=Dt*(sigq'*sigq + thetav)+(xtm1_tm1(end)-thetav)*(1-exp(-KAPPAv*Dt))/KAPPAv;
    OMEGA_xq=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))*SIGMA*sigq;
    OMEGA_v=thetav*SIGMAv^2*(1-exp(-KAPPAv*Dt))^2/(2*KAPPAv)+xtm1_tm1(end)*SIGMAv^2*(exp(-KAPPAv*Dt)-exp(-2*KAPPAv*Dt))/KAPPAv;
    OMEGA_vq=rho*SIGMAv*((1-exp(-KAPPAv/2*Dt))^2/KAPPAv*thetav+xtm1_tm1(end)*(exp(-KAPPAv/2*Dt)-exp(-KAPPAv*Dt))/(KAPPAv/2));
    OMEGA=[OMEGA_q, OMEGA_xq',OMEGA_vq;
           OMEGA_xq,OMEGA_x,zeros(Nfac,1);
           OMEGA_vq,zeros(1,Nfac),OMEGA_v];
       
    % prediction
    xt_tm1=A+B*xtm1_tm1;
    Pt_tm1=B*Ptm1_tm1*B'+OMEGA;

    % observation equation: y_t=a+b*z_t+eps^{y}_t
    if notmissing_p(j) == 0
        b= [zeros(Ny,1),by,cy];
        a= ay;
        y_obs=y_vec(j,:)';
        Rvec=delta_y.^2;
    elseif notmissing_p(j) == 1
        b= [1,zeros(1,Nfac),0;zeros(Ny,1),by,cy];
        a= [0;ay];
        y_obs=[p_vec(j); y_vec(j,:)'];
        Rvec=[delta_p^2; delta_y.^2];  
    end      

    if notmissing_bcf(j) == 1
        b = [b;[zeros(size(bf,1),1),bf,cf]];  
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

        b = [b;[0,bfLT,cfLT]];
        a = [a;afLT];

        y_obs= [y_obs; bcfLT_vec(j)];
        Rvec= [Rvec; delta_bcfLT^2];
    end  

    if notmissing_tips(j) == 1
        b = [b;[zeros(size(by_R,1),1),by_R,cy_R]];  
        a = [a;ay_R];
        y_obs= [y_obs; ytips_vec(j,:)'];
        Rvec= [Rvec; delta_tips.^2];
    end  


    if notmissing_caps(j)==1
        [ay_caps,by_caps,cy_caps,tmp_prc_caps] = CapsFacLoad_FRBAv_o(xt_tm1(2:Nfac+1),xt_tm1(Nfac+2), ...
            KAPPA,SIGMA,theta,rho0,rho1,rhov,lambda0,SIGMAlambda1,KAPPAv,SIGMAv,thetav, ...
            rho0_pi,rho1_pi,rhov_pi,sigq,gamv,gamvx,rho,MATgrid_caps,STRIKEgrid_caps);

        if isempty(ay_caps) | isempty(by_caps)
             warning('warning from caps pricing')
             logL=BadLike;
             logL_vec=[ ];
             state=[ ];
             model_prc_caps=[];
             model_prc_floors=[];
             return;
        end    
        model_prc_caps(j,:)=tmp_prc_caps';
        
        cy_caps=cy_caps./[1;1;1;3;3;3];
        by_caps=by_caps./repmat([1;1;1;3;3;3],1,Nfac);
        ay_caps=ay_caps./[1;1;1;3;3;3];
        caps_prc=caps_vec(j,:)'./[1;1;1;3;3;3];
        
        b = [b;[zeros(size(by_caps,1),1),by_caps,cy_caps]];  
        a = [a;ay_caps];
        y_obs = [y_obs; caps_prc];
        Rvec = [Rvec; delta_options^2*ones(Ncaps,1)];      
    end
    
    if notmissing_floors(j)==1
        [ay_floors,by_floors,cy_floors,tmp_prc_floors] = FloorsFacLoad_FRBAv_o(xt_tm1(2:Nfac+1),xt_tm1(Nfac+2), ...
            KAPPA,SIGMA,theta,rho0,rho1,rhov,lambda0,SIGMAlambda1,KAPPAv,SIGMAv,thetav, ...
            rho0_pi,rho1_pi,rhov_pi,sigq,gamv,gamvx,rho,MATgrid_floors,STRIKEgrid_floors);
                                
        if isempty(ay_floors) | isempty(by_floors) | isempty(cy_floors)
             warning('warning from floors pricing')
             logL=BadLike;
             logL_vec=[ ];
             state=[ ];
             model_prc_caps=[];
             model_prc_floors=[];
             return;
        end    
        model_prc_floors(j,:)=tmp_prc_floors';
        
        cy_floors=cy_floors./[1;1;1;3;3;3];
        by_floors=by_floors./repmat([1;1;1;3;3;3],1,Nfac);
        ay_floors=ay_floors./[1;1;1;3;3;3];
        floors_prc=floors_vec(j,:)'./[1;1;1;3;3;3];
        
        b = [b;[zeros(size(by_floors,1),1),by_floors,cy_floors]];  
        a = [a;ay_floors];
        y_obs = [y_obs; floors_prc];
        Rvec = [Rvec; delta_options^2*ones(Ncaps,1)];      
    end
    
    R=diag(Rvec);

    yt_tm1=a+b*xt_tm1;
    Vt_tm1=b*Pt_tm1*b'+R;    
    et=y_obs-yt_tm1;    % forecast error
    err(j,1:length(et))=et';

    rcVt=rcond(Vt_tm1);
    if (rcVt < 1.e-14 | isfinite(rcVt) == 0)
        warning('Vt is nearly singular')
        logL=9999;
        logL_vec=[ ];
        state=[ ];
        model_prc_caps=[];
        model_prc_floors=[];
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


