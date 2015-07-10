% function [logL,state,IE,IRP,IEx,IEv,IRPx,IRPv,model_caps_vec,model_floors_vec,FTP,FTP_R,IRPx_appr,IRPv_appr]=output_DKWoption_9013
% close all; clear; clc;

startdaten=datenum(1990,1,1);

% enddaten=datenum(2013,12,31);
% data_FRBA; 
% [Y,M,D]=datevec(mydate);
% load '.\results\pvar_DKWoption_9013.mat';
% id_lastobs=find(Y==2013 & M==12 & D==18);

enddaten=datenum(2012,12,31);
load '.\results\pvar_DKWoption_9012.mat';
data_FRBA; 
[Y,M,D]=datevec(mydate);
id_lastobs=find(Y==2012 & M==12 & D==19);



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

[logL,logL_vec,state,model_IE_options]= logL_DKWoption(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);

IE_options_avg=[mean(IE_options(:,1:3)')',mean(IE_options(:,4:6)')'];


figure; 
subplot(2,1,1);
plot([model_IE_options(id_options,1),IE_options_avg(id_options,1)])
subplot(2,1,2);
plot([model_IE_options(id_options,2),IE_options_avg(id_options,2)])



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


fig_state=figure; 
subplot(2,3,1); plot(Y+M/12+D/365,state(:,2),'LineWidth',4); title('DKW(options): state variable x1','FontSize',12);
subplot(2,3,2); plot(Y+M/12+D/365,state(:,3),'LineWidth',4); title('DKW(options): state variable x2','FontSize',12);
subplot(2,3,3); plot(Y+M/12+D/365,state(:,4),'LineWidth',4); title('DKW(options): state variable x3','FontSize',12);

% print(fig_state,'-depsc2','..\results\state_DKWoption.eps');



fig_IE=figure;
plot(Y+M/12+D/365,IE,'LineWidth',3); title('Inflation Expectations','FontSize',12); 
legend('1-year','10-year');

% print(fig_IE,'-depsc2','..\results\IE_DKWoption.eps');

load '..\data\umich.txt';
id_mich=find(umich(:,4)>0 & datenum(umich(:,1:3))>mydate(1) & datenum(umich(:,1:3))<=mydate(end));

load '..\data\FMAieS.txt';
id_ieS=find(FMAieS(:,4)>0 & datenum(FMAieS(:,1:3))>mydate(1) & datenum(FMAieS(:,1:3))<=mydate(end));


fig_allIE=figure;
plot(Y+M/12+D/365,IE(:,1),'b-'); hold on;
plot(umich(id_mich,1)+umich(id_mich,2)/12+umich(id_mich,3)/365,umich(id_mich,4)/100,'r-'); hold on;
plot(FMAieS(id_ieS,1)+FMAieS(id_ieS,2)/12+FMAieS(id_ieS,3)/365,FMAieS(id_ieS,4)/100,'ko'); hold off;
legend('model','Michigan','Blue Chip');
title('1-year Inflation Expectation');

% print(fig_allIE,'-depsc2','..\results\allIE_DKWoption.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample Inflation Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax=(eye(Nfac)-inv(KAPPA*Dt)*(eye(Nfac)-expm(-KAPPA*Dt)))*theta;
bx=(eye(Nfac)-expm(-KAPPA'*Dt))*inv(KAPPA'*Dt);

W=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))/Dt;


aI=rho0_pi+rho1_pi'*ax;
bI=rho1_pi'*bx; bI=bI';

hor_is=[1 6 12];    % forecasting horizon
Nhor_is=length(hor_is);

Np=length(id_p);
phat_is=NaN(Np,Nhor_is);
phat_is_RW=NaN(Np,Nhor_is);

p_is_rmse=zeros(Nhor_is,1);
p_is_RW_rmse=zeros(Nhor_is,1);

for j=1:Nhor_is
    p_prev=exp(p_vec(id_p(1)));
    hor=hor_is(j);
    
    for i=(hor+1):Np
        id_forecast=id_p(i-hor):(id_p(i)-1);

        phat_is(i,j)=p_prev*exp((aI*length(id_forecast)+sum(state(id_forecast,2:4))*bI)*Dt);

        phat_is_RW(i,j)=p_prev;

        p_prev=exp(p_vec(id_p(i-hor+1)));
    end
    fe_is=(phat_is((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
    fe_is_RW=(phat_is_RW((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
    
    p_is_rmse(j)=sqrt(mean((fe_is.^2)));
    p_is_RW_rmse(j)=sqrt(mean(fe_is_RW.^2));

end
[p_is_rmse,p_is_RW_rmse]

figure;
for j=1:3
    subplot(3,1,j);
    plot([exp(p_vec(id_p((1+hor_is(j)):end))),phat_is((1+hor_is(j)):end,j),phat_is_RW((1+hor_is(j)):end,j)])
    legend('data','model','rw');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample YIELD Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out-of-sample (OOS) Inflation Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startdaten=datenum(1990,1,1);
enddaten=datenum(2015,12,31);
data_FRBA; 
[Y,M,D]=datevec(mydate);
[logL,logL_vec,state,model_IE_options]= logL_DKWoption(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);

p_vec(1:(id_lastobs-1))=[];
state(1:(id_lastobs-1),:)=[];


hor_oos=[1 6 12];    % forecasting horizon in months
Nhor_oos=length(hor_oos);

id_p_oos=find(p_vec);
Np_oos=length(id_p_oos);

p_oos_rmse=zeros(Nhor_oos,1);
p_oos_RW_rmse=zeros(Nhor_oos,1);

for j=1:Nhor_oos
    p_prev=exp(p_vec(id_p_oos(1)));
    hor=hor_oos(j);

    phat_oos=NaN(Np_oos,1);
    phat_oos_RW=NaN(Np_oos,1);

    for i=(1+hor):Np_oos
        id_forecast=id_p_oos((i-hor)):(id_p_oos(i)-1);
        phat_oos(i)=p_prev*exp((aI*length(id_forecast)+sum(state(id_forecast,2:4))*bI)*Dt);

        phat_oos_RW(i)=p_prev;

        p_prev=exp(p_vec(id_p_oos(i-hor+1)));
    end

    fe_oos=(phat_oos((1+hor):Np_oos)-exp(p_vec(id_p_oos((1+hor):Np_oos))))./exp(p_vec(id_p_oos((1+hor):Np_oos)));
    fe_oos_RW=(phat_oos_RW((1+hor):Np_oos)-exp(p_vec(id_p_oos((1+hor):Np_oos))))./exp(p_vec(id_p_oos((1+hor):Np_oos)));

    p_oos_rmse(j)=sqrt(mean((fe_oos.^2)));
    p_oos_RW_rmse(j)=sqrt(mean(fe_oos_RW.^2));
end

[p_oos_rmse,p_oos_RW_rmse]

figure;
plot([exp(p_vec(id_p_oos)),phat_oos,phat_oos_RW])
legend('data','model','rw');




