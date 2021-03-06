function [logL,state,IE,IRP,model_caps_vec,model_floors_vec]=analysis_DKW_o
% close all; clear; clc;

is_year_beg=1990; % beginning year of in-sample period
is_year_end=2015; % ending year of in-sample period

startdaten=datenum(is_year_beg,1,1);
enddaten=datenum(is_year_end,12,31);
load '.\results\pvar_DKW_o_9012.mat';
data_FRBA; 
[Y,M,D]=datevec(mydate);

% model specification
KAPPA=NaN(Nfac);       
SIGMA=NaN(Nfac);
theta=NaN(Nfac,1);   
rho0=NaN;
rho1=NaN(Nfac,1);
lambda0=NaN(Nfac,1);
SIGMAlambda1=NaN(Nfac);
delta_y=NaN(Ny,1);
delta_bcf=NaN(2,1);
delta_bcfLT=NaN;

rho0_pi=NaN;
rho1_pi=NaN(Nfac,1);
rhov_pi=NaN;
sigq=NaN(Nfac,1);
sigqx=NaN;

delta_p=NaN;
delta_tips=NaN(Nytips,1);
delta_options=NaN(length(MATgrid_options),1);

% parameter normalization
KAPPA=diag(NaN(Nfac,1));
SIGMA(1,1)=0.01;
SIGMA(1,2)=0;
SIGMA(1,3)=0;
SIGMA(2,2)=0.01;
SIGMA(2,3)=0;
SIGMA(3,3)=0.01;

theta=zeros(Nfac,1);

delta_bcfLT=0.0075;

delta_p=0;

paras_vec=[reshape(KAPPA,Nfac^2,1);reshape(SIGMA,Nfac^2,1);theta;
      rho0;rho1;lambda0;reshape(SIGMAlambda1,Nfac^2,1);
      delta_y;delta_bcf;delta_bcfLT;rho0_pi;rho1_pi;sigq;
      sigqx;delta_p;delta_tips;delta_options];

j=1;
for i=1:length(paras_vec)
  if isnan(paras_vec(i))
    paras_vec(i)=pvar(j);
    pidx(j)=i;
    j=j+1;
  end
end
pidx=pidx'; 

                     
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

[logL,logL_vec,state,model_IE_options]= logL_DKW_o(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);


rho0_R=rho0-rho0_pi-1/2*(sigq'*sigq+sigqx^2)+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
lambda0_R=lambda0-sigq;
SIGMAlambda1_R=SIGMAlambda1;

%%%%%%%%%%%%%%%%%%%%%%%
% Yield Decomposition %
%%%%%%%%%%%%%%%%%%%%%%%

horLIST=[1 10];
IE=zeros(T,length(horLIST));
IRP=zeros(T,length(horLIST));

for i=1:length(horLIST);
    MAT=horLIST(i);
    [ay1,by1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MAT);

    [ay_R1,by_R1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0_R,rho1_R,lambda0_R,SIGMAlambda1_R,MAT);
    
    aI1=rho0_pi+rho1_pi'*(eye(Nfac)-1/MAT*inv(-KAPPA)*(expm(-KAPPA*MAT)-eye(Nfac)))*theta;
    bI1=1/MAT*inv(-KAPPA')*(expm(-KAPPA'*MAT)-eye(Nfac))*rho1_pi; 
    bI1=bI1';

    IE(:,i)=aI1+state(:,2:4)*bI1';
    IRP(:,i)=(ay1-ay_R1-aI1)+state(:,2:4)*(by1-by_R1-bI1)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-sample option prices fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATgrid_options=[1 3]; 
STRIKEgrid_options=[0.01 0.02 0.03];
model_caps_vec=zeros(T,length(MATgrid_options)*length(STRIKEgrid_options));
model_floors_vec=zeros(T,length(MATgrid_options)*length(STRIKEgrid_options));

for i=1:length(MATgrid_options)
    for j=1:length(STRIKEgrid_options)
        MAT=MATgrid_options(i);
        STRIKE=STRIKEgrid_options(j);
        tmpcount=j+length(STRIKEgrid_options)*(i-1);

        [ay1,by1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MAT);    
        tmp_Ay=-MAT*ay1; tmp_By=-MAT*by1';
        
        KAPPA_Q=KAPPA+SIGMAlambda1;
        theta_Q=inv(KAPPA_Q)*(KAPPA*theta-SIGMA*lambda0+SIGMA*SIGMA'*tmp_By);

        rho0_Q=rho0_pi-sigq'*lambda0+sigq'*SIGMA'*tmp_By;
        rho1_Q=rho1_pi-lambda1'*sigq;
                
        tmp_kron_Q=inv(kron(-KAPPA_Q,eye(Nfac))+kron(eye(Nfac),-KAPPA_Q));
        OMEGA_x_Q=reshape(tmp_kron_Q*reshape(expm(-KAPPA_Q*MAT)*SIGMA*SIGMA'*expm(-KAPPA_Q'*MAT)-SIGMA*SIGMA',Nfac^2,1),Nfac,Nfac);

        H0=((sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q)'*(sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q))*MAT ...
            -2*(sigq+SIGMA'*inv(KAPPA_Q')*rho1_Q)'*SIGMA'*inv(KAPPA_Q')*(eye(Nfac)-expm(-KAPPA_Q'*MAT))*inv(KAPPA_Q')*rho1_Q  ...
            +rho1_Q'*inv(KAPPA_Q)*OMEGA_x_Q*inv(KAPPA_Q')*rho1_Q;

        H1=MAT;
        
        tmpIE=(rho0_Q+rho1_Q'*theta_Q)+rho1_Q'*inv(KAPPA_Q*MAT) ...
              *(eye(Nfac)-expm(-KAPPA_Q*MAT))*(state(:,2:4)'-theta_Q*ones(1,T));
        tmpIE=tmpIE';
        tmpIU=1/MAT*(H0+sigqx^2*H1);
        
        tmph0=(-log(1+STRIKE)+tmpIE+tmpIU)./sqrt(tmpIU/MAT);
        tmph1=(-log(1+STRIKE)+tmpIE)./sqrt(tmpIU/MAT);
        
        
        model_caps_vec(:,tmpcount)= exp(tmp_Ay+state(:,2:4)*tmp_By).* ...
               (exp(MAT*(tmpIE+1/2*tmpIU)).*normcdf(tmph0)-(1+STRIKE)^MAT*normcdf(tmph1));
        model_floors_vec(:,tmpcount)= exp(tmp_Ay+state(:,2:4)*tmp_By).* ...
               (-exp(MAT*(tmpIE+1/2*tmpIU)).*normcdf(-tmph0)+(1+STRIKE)^MAT*normcdf(-tmph1));
    end
end

