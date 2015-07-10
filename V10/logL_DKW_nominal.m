function [logL,logL_vec,state]= logL_DKW_nominal(pvar,paras_vec,pidx,Nfac,dt, ...              
    y_vec,MATgrid,p_vec,notmissing_p)        

BadLike= 9999;

pall=pallGet(pvar,paras_vec,pidx);
[T,Ny]=size(y_vec);

KAPPA=reshape(pall(1:Nfac^2),Nfac,Nfac);
SIGMA=reshape(pall(1+Nfac^2:2*Nfac^2),Nfac,Nfac);
theta=pall(2*Nfac^2+1:2*Nfac^2+Nfac);
rho0=pall(2*Nfac^2+Nfac+1);
rho1=pall(2*Nfac^2+Nfac+2:2*Nfac^2+2*Nfac+1);
lambda0=pall(2*Nfac^2+2*Nfac+2:2*Nfac^2+3*Nfac+1);
SIGMAlambda1=reshape(pall(2*Nfac^2+3*Nfac+2:3*Nfac^2+3*Nfac+1),Nfac,Nfac);
delta_y=pall(3*Nfac^2+3*Nfac+2:3*Nfac^2+3*Nfac+1+Ny);

rho0_pi=pall(3*Nfac^2+3*Nfac+2+Ny);
rho1_pi=pall(3*Nfac^2+3*Nfac+3+Ny:3*Nfac^2+4*Nfac+2+Ny);
sigq=pall(3*Nfac^2+4*Nfac+3+Ny:3*Nfac^2+5*Nfac+2+Ny);
sigqx=pall(3*Nfac^2+5*Nfac+3+Ny);
delta_p=pall(3*Nfac^2+5*Nfac+4+Ny);

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
  
% x_{t+1}=Ax+Bx*x_t+eps^{x}_{t+1}
Bx=expm(-KAPPA*dt);
Ax= (eye(Nfac)-Bx)*theta;

B=zeros(Nfac+1);
B(1,1)=1; B(1,2:end)=rho1_pi'*dt;
B(2:end,2:end)=Bx;

A=zeros(Nfac+1,1);
A(1)=rho0_pi*dt;
A(2:end)=Ax;

tmp_kron=inv(kron(-KAPPA,eye(Nfac))+kron(eye(Nfac),-KAPPA));
OMEGA_x_ss=-reshape(tmp_kron*reshape(SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);  
OMEGA_x=reshape(tmp_kron*reshape(expm(-KAPPA*dt)*SIGMA*SIGMA'*expm(-KAPPA'*dt)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

OMEGA_ss=zeros(Nfac+1);
OMEGA_ss(1,1)=1e6;  % diffuse prior for price: choose a large enough number
OMEGA_ss(2:end,2:end)=OMEGA_x_ss;

% OMEGA: var[x_t|(t-1)]
OMEGA_q=dt*(sigq'*sigq + sigqx^2);
OMEGA_xq=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*dt))*SIGMA*sigq;
OMEGA=[OMEGA_q, OMEGA_xq';
       OMEGA_xq,OMEGA_x];

ik=find(p_vec(1:30));
priceinit=p_vec(ik(1));

% unconditional mean and variance
xtm1_tm1=[priceinit;theta];         % E(xt|t)
Ptm1_tm1=OMEGA_ss;                  % V(xt|t)


[T,Ny]=size(y_vec);
logL_vec=zeros(T,1);
state=zeros(T,Nfac+1);
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
   
    R=diag(Rvec);

    yt_tm1=a+b*xt_tm1;
    Vt_tm1=b*Pt_tm1*b'+R;    
    et=y_obs-yt_tm1;    % forecast error
  
    rcVt=rcond(Vt_tm1);
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


