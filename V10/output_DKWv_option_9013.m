% function [logL,state,IE,IRP,IEx,IEv,IRPx,IRPv,model_caps_vec,model_floors_vec,FTP,FTP_R,IRPx_appr,IRPv_appr]=output_DKWv_option_9013
% close all; clear; clc;

% startdaten=datenum(1990,1,1);
% enddaten=datenum(2013,12,31);
% data_FRBA; 
% [Y,M,D]=datevec(mydate);
% load '.\results\pvar_DKWoption_9013.mat';
% id_lastobs=find(Y==2013 & M==12 & D==18);

is_year_beg=1990; % beginning year of in-sample period
is_year_end=2012; % ending year of in-sample period

startdaten=datenum(is_year_beg,1,1);
enddaten=datenum(is_year_end,12,31);
load '.\results\pvar_DKWv_option_9012.mat';
data_FRBA; 
[Y,M,D]=datevec(mydate);
id_lastobs=find(Y==2012 & M==12 & D==19);



% model specification
KAPPA=NaN(Nfac);       
SIGMA=NaN(Nfac);
theta=NaN(Nfac,1);   
rho0=NaN;
rho1=NaN(Nfac,1);
rhov=NaN;
lambda0=NaN(Nfac,1);
SIGMAlambda1=NaN(Nfac);
delta_y=NaN(Ny,1);
delta_bcf=NaN(2,1);
delta_bcfLT=NaN;

rho0_pi=NaN;
rho1_pi=NaN(Nfac,1);
rhov_pi=NaN;
sigq=NaN(Nfac,1);

KAPPAv=NaN;
SIGMAv=NaN;
thetav=NaN;
gamv=NaN;
gamvx=NaN;
rho=NaN;

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
      rho0;rho1;rhov;lambda0;reshape(SIGMAlambda1,Nfac^2,1);
      delta_y;delta_bcf;delta_bcfLT;rho0_pi;rho1_pi;rhov_pi;sigq;
      KAPPAv;SIGMAv;thetav;gamv;gamvx;rho;delta_p;delta_tips;delta_options];

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


lambda1=inv(SIGMA)*SIGMAlambda1;

[logL,logL_vec,state,model_IE_options]= logL_DKWv_option(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);

IE_options_avg=[mean(IE_options(:,1:3)')',mean(IE_options(:,4:6)')'];


figure; 
subplot(2,1,1);
plot([model_IE_options(id_options,1),IE_options_avg(id_options,1)])
subplot(2,1,2);
plot([model_IE_options(id_options,2),IE_options_avg(id_options,2)])



rho0_R=rho0-rho0_pi-1/2*sigq'*sigq+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
rhov_R=rhov-rhov_pi+gamvx-1/2;
lambda0_R=lambda0-sigq;
lambda1_R=lambda1;
SIGMAlambda1_R=SIGMA*lambda1_R;
KAPPA_rn_R=KAPPA + SIGMAlambda1_R;
KAPPAtheta_rn_R=KAPPA*theta - SIGMA*lambda0_R;

%%%%%%%%%%%%%%%%%%%%%%%
% Yield Decomposition %
%%%%%%%%%%%%%%%%%%%%%%%

horLIST=[1 10];
IE=zeros(T,length(horLIST));
IRP=zeros(T,length(horLIST));

for i=1:length(horLIST);
    MAT=horLIST(i);
    ax=(eye(Nfac)-inv(KAPPA*MAT)*(eye(Nfac)-expm(-KAPPA*MAT)))*theta;
    bx=(eye(Nfac)-expm(-KAPPA'*MAT))*inv(KAPPA'*MAT);
    av=(1-inv(KAPPAv*MAT)*(1-exp(-KAPPAv*MAT)))*thetav;
    bv=(1-exp(-KAPPAv'*MAT))*inv(KAPPAv'*MAT);
    
    aI1=rho0_pi+rho1_pi'*ax+rhov_pi*av;
    bI1=bx'*rho1_pi; bI1=bI1';
    cI1=rhov_pi*bv;    

    KAPPA_rn=KAPPA + SIGMAlambda1;
    KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;
    KAPPAv_rn=KAPPAv+sqrt(1-rho^2)*SIGMAv*gamv+rho*SIGMAv*gamvx;    
    
    [ay,by,cy] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA*SIGMA',rho0,rho1,rhov,KAPPAv_rn,SIGMAv^2);

    ay=-ay./MAT; ay=ay';
    by=-by./repmat(MAT,Nfac,1); by=by';
    cy=-cy./MAT; cy=cy';
    ay1=ay; by1=by; cy1=cy;

    rho0_R=rho0-rho0_pi-1/2*sigq'*sigq+lambda0'*sigq;
    rho1_R=rho1-rho1_pi+lambda1'*sigq;
    rhov_R=rhov-rhov_pi+gamvx-1/2;
    lambda0_R=lambda0-sigq;
    lambda1_R=lambda1;
    SIGMAlambda1_R=SIGMA*lambda1_R;
    KAPPA_rn_R=KAPPA + SIGMAlambda1_R;
    KAPPAtheta_rn_R=KAPPA*theta - SIGMA*lambda0_R;
    [ay_R,by_R,cy_R] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA',rho0_R,rho1_R,rhov_R,KAPPAv_rn,SIGMAv^2);
    ay_R=-ay_R./MAT; ay_R=ay_R';
    by_R=-by_R./repmat(MAT,Nfac,1); by_R=by_R';
    cy_R=-cy_R./MAT; cy_R=cy_R';
    ay_R1=ay_R; by_R1=by_R; cy_R1=cy_R;


    IE(:,i)=aI1+state(:,2:4)*bI1'+state(:,5)*cI1;
    IRP(:,i)=(ay1-ay_R1-aI1)+state(:,2:4)*(by1-by_R1-bI1)'+state(:,5)*(cy1-cy_R1-cI1)';

end


fig_state=figure; 
subplot(2,3,1); plot(Y+M/12+D/365,state(:,2),'LineWidth',4); title('DKWv(options): state variable x1','FontSize',12);
subplot(2,3,2); plot(Y+M/12+D/365,state(:,3),'LineWidth',4); title('DKWv(options): state variable x2','FontSize',12);
subplot(2,3,3); plot(Y+M/12+D/365,state(:,4),'LineWidth',4); title('DKWv(options): state variable x3','FontSize',12);
subplot(2,3,4); plot(Y+M/12+D/365,state(:,5),'LineWidth',4); title('DKWv(options): state variable v','FontSize',12);

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
% Various In-Sample Forecasts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample Nominal Yields Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT_allyields=[0.25 0.5 1:30];
MAT_is=[0.5 1 2 3 4 5]; % bond maturity in years (6m/2y/10y)
NMAT_is=length(MAT_is);
hor_is=[13 26 52]; % forecasting horizon in weeks (3m/6m/12m)
Nhor_is=length(hor_is);

yN_is_rmse=zeros(NMAT_is,Nhor_is);
yN_is_RW_rmse=zeros(NMAT_is,Nhor_is);

for k=1:NMAT_is
    MAT=MAT_is(k);
    id_allyields=find(MAT_allyields==MAT);
        
    [ay,by,cy] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA*SIGMA',rho0,rho1,rhov,KAPPAv_rn,SIGMAv^2);

    ay=-ay./MAT; ay=ay';
    by=-by./repmat(MAT,Nfac,1); by=by';
    cy=-cy./MAT; cy=cy';
    ay1=ay; by1=by; cy1=cy;
    
    for j=1:Nhor_is
        hor=hor_is(j);

        yNhat_is=NaN(T,1);
        yNhat_is_RW=NaN(T,1);

        for i=(hor+1):T
            yNhat_is(i)=ay1+by1*(theta+expm(-KAPPA*hor*Dt)*(state(i-hor,2:4)'-theta)) ...
                           +cy1*(thetav+exp(-KAPPAv*hor*Dt)*(state(i-hor,5)-thetav));
            yNhat_is_RW(i)=yields(i-hor,id_allyields);
        end
        
        fe_is=(yNhat_is((hor+1):T)-yields((hor+1):T,id_allyields));
        fe_is_RW=(yNhat_is_RW((hor+1):T)-yields((hor+1):T,id_allyields));        
        
        yN_is_rmse(k,j)=sqrt(mean((fe_is.^2)));
        yN_is_RW_rmse(k,j)=sqrt(mean(fe_is_RW.^2));        
    end    
end
MAT_list=repmat(MAT_is',Nhor_is,1);
HOR_list=repmat(floor(hor_is/4.3),NMAT_is,1); HOR_list=HOR_list(:);
[MAT_list,HOR_list,yN_is_rmse(:)*100,yN_is_RW_rmse(:)*100]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample Real Yields Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT_alltips=[2:30];
MAT_is=[2 3 4 5]; % bond maturity in years (6m/2y/10y)
NMAT_is=length(MAT_is);
hor_is=[13 26 52]; % forecasting horizon in weeks (3m/6m/12m)
Nhor_is=length(hor_is);

yR_is_rmse=zeros(NMAT_is,Nhor_is);
yR_is_RW_rmse=zeros(NMAT_is,Nhor_is);

id_alltips_start=min(find(alltips(:,1)));

for k=1:NMAT_is
    MAT=MAT_is(k);
    id_alltips=find(MAT_alltips==MAT);
        
    [ay_R,by_R,cy_R] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA',rho0_R,rho1_R,rhov_R,KAPPAv_rn,SIGMAv^2);
    ay_R=-ay_R./MAT; ay_R=ay_R';
    by_R=-by_R./repmat(MAT,Nfac,1); by_R=by_R';
    cy_R=-cy_R./MAT; cy_R=cy_R';
    ay_R1=ay_R; by_R1=by_R; cy_R1=cy_R;
    
    for j=1:Nhor_is
        hor=hor_is(j);

        yRhat_is=NaN(T,1);
        yRhat_is_RW=NaN(T,1);

        for i=(hor+id_alltips_start):T
            yRhat_is(i)=ay_R1+by_R1*(theta+expm(-KAPPA*hor*Dt)*(state(i-hor,2:4)'-theta)) ...
                           +cy_R1*(thetav+exp(-KAPPAv*hor*Dt)*(state(i-hor,5)-thetav));
            yRhat_is_RW(i)=alltips(i-hor,id_alltips);
        end
        
        fe_is=(yRhat_is((hor+id_alltips_start):T)-alltips((hor+id_alltips_start):T,id_alltips));
        fe_is_RW=(yRhat_is_RW((hor+id_alltips_start):T)-alltips((hor+id_alltips_start):T,id_alltips));        
        
        yR_is_rmse(k,j)=sqrt(mean((fe_is.^2)));
        yR_is_RW_rmse(k,j)=sqrt(mean(fe_is_RW.^2));        
    end    
end
MAT_list=repmat(MAT_is',Nhor_is,1);
HOR_list=repmat(floor(hor_is/4.3),NMAT_is,1); HOR_list=HOR_list(:);
[MAT_list,HOR_list,yR_is_rmse(:)*100,yR_is_RW_rmse(:)*100]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample Inflation Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax=(eye(Nfac)-inv(KAPPA*Dt)*(eye(Nfac)-expm(-KAPPA*Dt)))*theta;
bx=(eye(Nfac)-expm(-KAPPA'*Dt))*inv(KAPPA'*Dt);

av=(1-inv(KAPPAv*Dt)*(1-exp(-KAPPAv*Dt)))*thetav;
bv=(1-exp(-KAPPAv'*Dt))*inv(KAPPAv'*Dt);

aI=rho0_pi+rho1_pi'*ax+rhov_pi*av;
bI=rho1_pi'*bx; bI=bI';
cI=rhov_pi*bv;


hor_is=[1 6 12];    % forecasting horizon (months)
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

        phat_is(i,j)=p_prev*exp((aI*length(id_forecast)+sum(state(id_forecast,2:4))*bI+sum(state(id_forecast,5))*cI)*Dt);

        phat_is_RW(i,j)=p_prev;

        p_prev=exp(p_vec(id_p(i-hor+1)));
    end
    fe_is=(phat_is((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
    fe_is_RW=(phat_is_RW((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
    
    p_is_rmse(j)=sqrt(mean((fe_is.^2)));
    p_is_RW_rmse(j)=sqrt(mean(fe_is_RW.^2));

end
[hor_is',p_is_rmse,p_is_RW_rmse]

figure;
for j=1:3
    subplot(3,1,j);
    plot([exp(p_vec(id_p((1+hor_is(j)):end))),phat_is((1+hor_is(j)):end,j),phat_is_RW((1+hor_is(j)):end,j)])
    legend('data','model','rw');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Various Out-of-Sample Forecasts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startdaten=datenum(1990,1,1);
enddaten=datenum(2015,12,31);
data_FRBA; 
[Y,M,D]=datevec(mydate);
[logL,logL_vec,state,model_IE_options]= logL_DKWv_option(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out-of-Sample Nominal Yields Forecasting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_outsample=find(Y>is_year_end);
state_oos=state(id_outsample,:);

MAT_oos=[0.5 1 2 3 4 5]; % bond maturity in years (6m/2y/10y)
NMAT_oos=length(MAT_oos);
hor_oos=[13 26 52]; % forecasting horizon in weeks (3m/6m/12m)
Nhor_oos=length(hor_oos);

T_oos=length(id_outsample);
yields_oos=yields(id_outsample,:);

yN_oos_rmse=zeros(NMAT_oos,Nhor_oos);
yN_oos_RW_rmse=zeros(NMAT_oos,Nhor_oos);

for k=1:NMAT_oos
    MAT=MAT_oos(k);
    id_allyields=find(MATgrid_allyield==MAT);
        
    [ay,by,cy] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn,KAPPA_rn,KAPPAv*thetav,SIGMA*SIGMA',rho0,rho1,rhov,KAPPAv_rn,SIGMAv^2);

    ay=-ay./MAT; ay=ay';
    by=-by./repmat(MAT,Nfac,1); by=by';
    cy=-cy./MAT; cy=cy';
    ay1=ay; by1=by; cy1=cy;

        
    for j=1:Nhor_oos
        hor=hor_oos(j);
    
        yNhat_oos=NaN(T_oos,1);
        yNhat_oos_RW=NaN(T_oos,1);

        for i=(hor+1):T_oos
            yNhat_oos(i)=ay1+by1*(theta+expm(-KAPPA*hor*Dt)*(state_oos(i-hor,2:4)'-theta)) ...
                            +cy1*(thetav+exp(-KAPPAv*hor*Dt)*(state_oos(i-hor,5)-thetav));
            yNhat_oos_RW(i)=yields_oos(i-hor,id_allyields);
        end
        
        fe_oos=(yNhat_oos((hor+1):T_oos)-yields_oos((hor+1):T_oos,id_allyields));
        fe_oos_RW=(yNhat_oos_RW((hor+1):T_oos)-yields_oos((hor+1):T_oos,id_allyields));        
        
        yN_oos_rmse(k,j)=sqrt(mean((fe_oos.^2)));
        yN_oos_RW_rmse(k,j)=sqrt(mean(fe_oos_RW.^2));        
    end    
end
MAT_list=repmat(MAT_oos',Nhor_oos,1);
HOR_list=repmat(floor(hor_oos/4.3),NMAT_oos,1); HOR_list=HOR_list(:);
[MAT_list,HOR_list,yN_oos_rmse(:)*100,yN_oos_RW_rmse(:)*100]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out-of-Sample REAL Yields Forecasting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT_oos=[2 3 4 5]; % bond maturity in years (6m/2y/10y)
NMAT_oos=length(MAT_oos);
hor_oos=[13 26 52]; % forecasting horizon in weeks (3m/6m/12m)
Nhor_oos=length(hor_oos);

T_oos=length(id_outsample);
alltips_oos=alltips(id_outsample,:);

yR_oos_rmse=zeros(NMAT_oos,Nhor_oos);
yR_oos_RW_rmse=zeros(NMAT_oos,Nhor_oos);

for k=1:NMAT_oos
    MAT=MAT_oos(k);
    id_alltips=find(MATgrid_alltips==MAT);
        
    [ay_R,by_R,cy_R] = YieldFacLoad_ODE3(MAT,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA',rho0_R,rho1_R,rhov_R,KAPPAv_rn,SIGMAv^2);
    ay_R=-ay_R./MAT; ay_R=ay_R';
    by_R=-by_R./repmat(MAT,Nfac,1); by_R=by_R';
    cy_R=-cy_R./MAT; cy_R=cy_R';
    ay_R1=ay_R; by_R1=by_R; cy_R1=cy_R;
        
    for j=1:Nhor_oos
        hor=hor_oos(j);
    
        yRhat_oos=NaN(T_oos,1);
        yRhat_oos_RW=NaN(T_oos,1);

        for i=(hor+1):T_oos
            yRhat_oos(i)=ay_R1+by_R1*(theta+expm(-KAPPA*hor*Dt)*(state_oos(i-hor,2:4)'-theta)) ...
                              +cy_R1*(thetav+exp(-KAPPAv*hor*Dt)*(state_oos(i-hor,5)-thetav));
            yRhat_oos_RW(i)=alltips_oos(i-hor,id_alltips);
        end
        
        fe_oos=(yRhat_oos((hor+1):T_oos)-alltips_oos((hor+1):T_oos,id_allyields));
        fe_oos_RW=(yRhat_oos_RW((hor+1):T_oos)-alltips_oos((hor+1):T_oos,id_allyields));        
        
        yR_oos_rmse(k,j)=sqrt(mean((fe_oos.^2)));
        yR_oos_RW_rmse(k,j)=sqrt(mean(fe_oos_RW.^2));        
    end    
end
MAT_list=repmat(MAT_oos',Nhor_oos,1);
HOR_list=repmat(floor(hor_oos/4.3),NMAT_oos,1); HOR_list=HOR_list(:);
[MAT_list,HOR_list,yR_oos_rmse(:)*100,yR_oos_RW_rmse(:)*100]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out-of-sample (OOS) Inflation Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        phat_oos(i)=p_prev*exp((aI*length(id_forecast)+sum(state(id_forecast,2:4))*bI+sum(state(id_forecast,5))*cI)*Dt);

        phat_oos_RW(i)=p_prev;

        p_prev=exp(p_vec(id_p_oos(i-hor+1)));
    end

    fe_oos=(phat_oos((1+hor):Np_oos)-exp(p_vec(id_p_oos((1+hor):Np_oos))))./exp(p_vec(id_p_oos((1+hor):Np_oos)));
    fe_oos_RW=(phat_oos_RW((1+hor):Np_oos)-exp(p_vec(id_p_oos((1+hor):Np_oos))))./exp(p_vec(id_p_oos((1+hor):Np_oos)));

    p_oos_rmse(j)=sqrt(mean((fe_oos.^2)));
    p_oos_RW_rmse(j)=sqrt(mean(fe_oos_RW.^2));
end

[hor_oos',p_oos_rmse,p_oos_RW_rmse]

figure;
plot([exp(p_vec(id_p_oos)),phat_oos,phat_oos_RW])
legend('data','model','rw');




