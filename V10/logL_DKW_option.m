function [logL,logL_vec,state,model_IE_options]= logL_DKW_option(pvar,paras_vec,pidx,Nfac,Dt, ...              
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
delta_options=pall(3*Nfac^2+5*Nfac+8+Ny+Nytips:3*Nfac^2+5*Nfac+7+Ny+Nytips+length(MATgrid_options));

lambda1=inv(SIGMA)*SIGMAlambda1;

if min(diag(KAPPA)) < 2.e-4  
    warning('KAPPA is too singular')
    logL=9999;
    logL_vec=[];
    state=[];
    return;
end
% Compute the factor loadings
[ay,by]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MATgrid);
  
rho0_R=rho0-rho0_pi-1/2*(sigq'*sigq+sigqx^2)+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
lambda0_R=lambda0-sigq;
SIGMAlambda1_R=SIGMAlambda1;

[ay_R,by_R]=YieldFacLoad(KAPPA,SIGMA,theta,rho0_R,rho1_R,lambda0_R,SIGMAlambda1_R,TIPSgrid);

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

B=zeros(Nfac+1);
B(1,1)=1; B(1,2:end)=rho1_pi'*Dt;
B(2:end,2:end)=Bx;

A=zeros(Nfac+1,1);
A(1)=rho0_pi*Dt;
A(2:end)=Ax;

tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x_ss=-reshape(tmp_kron*reshape(SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);  
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*Dt)*SIGMA*SIGMA'*expm(-KAPPA'*Dt)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

aI_Q=zeros(length(MATgrid_options),1);
bI_Q=zeros(length(MATgrid_options),Nfac);

for i=1:length(MATgrid_options)
    MAT=MATgrid_options(i);
%     [aI,bI] = InfExpFacLoad(KAPPA,SIGMA,theta,sigq,sigqx,rho0_pi,rho1_pi,MAT);
    
    [tmp_ay,tmp_by]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MAT);

    KAPPA_Q=KAPPA+SIGMAlambda1;
    theta_Q=inv(KAPPA_Q)*(KAPPA*theta-SIGMA*lambda0+SIGMA*SIGMA'*tmp_by');
    
    rho0_Q=rho0_pi-lambda0'*sigq+sigq'*SIGMA'*tmp_by';
    rho1_Q=rho1_pi-lambda1'*sigq;
        
    [tmp_aI_Q,tmp_bI_Q] = InfExpFacLoad(KAPPA_Q,SIGMA,theta_Q,sigq,sigqx,rho0_Q,rho1_Q,MAT);

    aI_Q(i)=tmp_aI_Q;
    bI_Q(i,:)=tmp_bI_Q;    
end

OMEGA_ss=zeros(Nfac+1);
OMEGA_ss(1,1)=1e6;  % diffuse prior for price: choose a large enough number
OMEGA_ss(2:end,2:end)=OMEGA_x_ss;

% OMEGA: var[x_t|(t-1)]
OMEGA_q=Dt*(sigq'*sigq + sigqx^2);
OMEGA_xq=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))*SIGMA*sigq;
OMEGA=[OMEGA_q, OMEGA_xq';
       OMEGA_xq,OMEGA_x];


IE_options_avg=[mean(IE_options(:,1:3)')',mean(IE_options(:,4:6)')'];

   
ik=find(p_vec(1:30));
priceinit=p_vec(ik(1));

% unconditional mean and variance
xtm1_tm1=[priceinit;theta];         % E(xt|t)
Ptm1_tm1=OMEGA_ss;                  % V(xt|t)


[T,Ny]=size(y_vec);
logL_vec=zeros(T,1);
state=zeros(T,Nfac+1);
model_IE_options=zeros(T,length(MATgrid_options));

err=zeros(T,Ny+Nytips+Nhor+2);
for j=1:T  
    % prediction
    xt_tm1=A+B*xtm1_tm1;
    Pt_tm1=B*Ptm1_tm1*B'+OMEGA;

    % observation equation: y_t=a+b*x_t+eps^{y}_t
    if notmissing_p(j) == 0
        b= [zeros(Ny,1),by];
        a= ay;
        y_obs=y_vec(j,:)';
        Rvec=delta_y.^2;
    elseif notmissing_p(j) == 1
        b= [1,zeros(1,Nfac);zeros(Ny,1),by];
        a= [0;ay];
        y_obs=[p_vec(j); y_vec(j,:)'];
        Rvec=[delta_p^2; delta_y.^2];  
    end      
  
    if notmissing_bcf(j) == 1
        b = [b;[zeros(size(bf,1),1),bf]];  
        a = [a;af];
        y_obs= [y_obs; bcf_vec(j,:)'];
        Rvec= [Rvec; delta_bcf.^2];
    end  

    if notmissing_bcfLT(j) == 1      
        hor=horLT(j);
        W=inv(KAPPA)*(expm(-KAPPA*hor)-expm(-KAPPA*(hor+5.0)))/5.0;
        afLT= ay3 + by3*(eye(Nfac)-W)*theta;
        bfLT= by3*W;

        b = [b;[0,bfLT]];
        a = [a;afLT];

        y_obs= [y_obs; bcfLT_vec(j)];
        Rvec= [Rvec; delta_bcfLT^2];
    end  

    if notmissing_tips(j) == 1
        b = [b;[zeros(size(by_R,1),1),by_R]];  
        a = [a;ay_R];
        y_obs= [y_obs; ytips_vec(j,:)'];
        Rvec= [Rvec; delta_tips.^2];
    end  
    
    if notmissing_options(j) == 1
        b = [b;[zeros(size(bI_Q,1),1),bI_Q]];  
        a = [a;aI_Q];
        y_obs= [y_obs; IE_options_avg(j,:)'];
        Rvec= [Rvec; delta_options.^2];
        
        model_IE_options(j,:)=(aI_Q+bI_Q*xt_tm1(2:4))';
    end  
    
    R=diag(Rvec);

    yt_tm1=a+b*xt_tm1;
    Vt_tm1=b*Pt_tm1*b'+R;    
    et=y_obs-yt_tm1;    % forecast error
    err(j,1:length(et))=et';

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


