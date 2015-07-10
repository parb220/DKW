function [logL,state,IE,IRP,IEx,IEv,IRPx,IRPv,model_caps_vec,model_floors_vec,FTP,FTP_R,IRPx_appr,IRPv_appr]=output_FRBAoil_9013
% close all; clear; clc;

startdaten=datenum(1990,1,1);
% startdaten=datenum(2010,2,17);
enddaten=datenum(2013,12,31);

load '.\results\pvar_FRBAoil_9013.mat';

data_FRBA; 
[Y,M,D]=datevec(mydate);

% model specification
KAPPA=NaN(Nfac);       
SIGMA=NaN(Nfac);
theta=NaN(Nfac,1);   
rho0=NaN;
rho1=NaN(Nfac,1);
rhod=NaN;
rhos=NaN;
lambda0=NaN(Nfac,1);
SIGMAlambda1=NaN(Nfac);
delta_y=NaN(Ny,1);
delta_bcf=NaN(2,1);
delta_bcfLT=NaN;

rho0_pi=NaN;
rho1_pi=NaN(Nfac,1);
rhod_pi=NaN;
rhos_pi=NaN;

sigq=NaN(Nfac,1);
sigqx=NaN;

SIGMAs=NaN;
KAPPAd=NaN;
SIGMAd=NaN;
thetad=NaN;
lambda0_d=NaN;
lambda1_d=NaN;

delta_p=NaN;
delta_tips=NaN(Nytips,1);
delta_oil=NaN;

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
      rho0;rho1;rhod;rhos;lambda0;reshape(SIGMAlambda1,Nfac^2,1);
      delta_y;delta_bcf;delta_bcfLT;rho0_pi;rho1_pi;rhod_pi;rhos_pi;sigq;sigqx;
      SIGMAs;KAPPAd;SIGMAd;thetad;lambda0_d;lambda1_d;
      delta_p;delta_tips;delta_oil];

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


[logL,logL_vec,state,model_oil_vec]= logL_FRBAoil(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,oil_vec,MATgrid_oil,notmissing_oil);



%%%%%%%%%%%%%%%%%%%%%%%
% Yield Decomposition %
%%%%%%%%%%%%%%%%%%%%%%%

horLIST=[1 10];
IE=zeros(T,length(horLIST));
IRP=zeros(T,length(horLIST));

IEx=zeros(T,length(horLIST));
IEv=zeros(T,length(horLIST));

% IRPx_appr=zeros(T,length(horLIST));
% IRPv_appr=zeros(T,length(horLIST));

IRPx=zeros(T,length(horLIST));
IRPv=zeros(T,length(horLIST));

for i=1:length(horLIST);
    MAT=horLIST(i);
    
    ax=(eye(Nfac)-inv(KAPPA*MAT)*(eye(Nfac)-expm(-KAPPA*MAT)))*theta;
    bx=(eye(Nfac)-expm(-KAPPA'*MAT))*inv(KAPPA'*MAT);
    
    ad=(1-inv(KAPPAd*MAT)*(1-exp(-KAPPAd*MAT)))*thetad;
    bd=(1-exp(-KAPPAd'*MAT))*inv(KAPPAd'*MAT);
        
    W=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*MAT))/MAT;
    Wd=inv(KAPPAd)*(1-exp(-KAPPAd*MAT))/MAT;
    Ws=inv(phis)*(exp(phis*MAT)-1)/MAT;

    bsx=phi1'*inv(KAPPA+phis*eye(Nfac))*(Ws*eye(Nfac)-W);
    bsd=phid*inv(KAPPAd+phis)*(Ws-Wd);        
    bss=Ws;
    as=1/phis*(Ws-1)*(phi0+phi1'*theta+phid*thetad) ...
         -bsx*theta-bsd*thetad; 
        
    aI1=rho0_pi+rho1_pi'*ax+rhod_pi*ad+rhos_pi*as;
    bI1=rho1_pi'*bx+rhos_pi*bsx; bI1=bI1';
    cI2=rhod_pi*bd+rhos_pi*bsd;
    cI3=rhos_pi*bss;
    

    KAPPA_rn=KAPPA + SIGMAlambda1;
    KAPPAtheta_rn=KAPPA*theta - SIGMA*lambda0;
    theta_rn=inv(KAPPA_rn)*KAPPAtheta_rn;

    KAPPAd_rn=KAPPAd+SIGMAd*lambda1_d;
    thetad_rn=(KAPPAd*thetad-SIGMAd*lambda0_d)/KAPPAd_rn;

    phi0_rn=rho0-1/2*SIGMAs^2;
    phi1_rn=rho1;
    phid_rn=rhod-1;
    phis_rn=rhos;


    % Compute the nominal yield factor loadings
    tmp_k=[KAPPAtheta_rn;KAPPAd_rn*thetad_rn;phi0_rn];
    tmp_K=[-KAPPA_rn',zeros(Nfac,1),phi1_rn;zeros(1,Nfac),-KAPPAd_rn,phid_rn;zeros(1,Nfac),0,phis_rn];
    tmp_K=(-tmp_K)';
    tmp_H=zeros(Nfac+2,Nfac+2); 
    tmp_H(1:Nfac,1:Nfac)=SIGMA*SIGMA'; tmp_H(Nfac+1,Nfac+1)=SIGMAd^2; tmp_H(Nfac+2,Nfac+2)=SIGMAs^2;
    [ay_tran,by_tran] = YieldFacLoad_ODE(MAT,tmp_k,tmp_K,tmp_H,rho0,[rho1;rhod;rhos]);

    ay=ay_tran; by=by_tran(1:Nfac,:); dy=by_tran(Nfac+1,:); ey=by_tran(Nfac+2,:);

    ay=-ay./MAT; ay=ay';
    by=-by./repmat(MAT,Nfac,1); by=by';
    dy=-dy./MAT; dy=dy';
    ey=-ey./MAT; ey=ey';


    ay1=ay; by1=by; dy1=dy; ey1=ey;

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
    [ay_tran_R,by_tran_R] = YieldFacLoad_ODE(MAT,tmp_k_R,tmp_K_R,tmp_H,rho0_R,[rho1_R;rhod_R;rhos_R]);
    ay_R=ay_tran_R; by_R=by_tran_R(1:Nfac,:); dy_R=by_tran_R(Nfac+1,:); ey_R=by_tran_R(Nfac+2,:);

    ay_R=-ay_R./MAT; ay_R=ay_R';
    by_R=-by_R./repmat(MAT,Nfac,1); by_R=by_R';
    dy_R=-dy_R./MAT; dy_R=dy_R';
    ey_R=-ey_R./MAT; ey_R=ey_R';

    ay_R1=ay_R; by_R1=by_R; dy_R1=dy_R; ey_R1=ey_R;

    IE(:,i)=aI1+state(:,2:4)*bI1+state(:,5)*cI2+state(:,6)*cI3;
    IRP(:,i)=(ay1-ay_R1-aI1)+state(:,2:4)*(by1-by_R1-bI1')'+state(:,5)*(dy1-dy_R1-cI2)'+state(:,6)*(ey1-ey_R1-cI3)';

end


fig_state=figure; 
subplot(2,3,1); plot(Y+M/12+D/365,state(:,2),'LineWidth',4); title('FRBA (flex caps): state variable x1','FontSize',12);
subplot(2,3,2); plot(Y+M/12+D/365,state(:,3),'LineWidth',4); title('FRBA (flex caps): state variable x2','FontSize',12);
subplot(2,3,3); plot(Y+M/12+D/365,state(:,4),'LineWidth',4); title('FRBA (flex caps): state variable x3','FontSize',12);
subplot(2,3,4); plot(Y+M/12+D/365,state(:,5),'LineWidth',4); title('FRBA (flex caps): state variable delta','FontSize',12);  
subplot(2,3,5); plot(Y+M/12+D/365,state(:,6),'LineWidth',4); title('FRBA (flex caps): state variable spot price','FontSize',12); 

% print(fig_state,'-depsc2','..\results\state_FRBAoil.eps');



id_oil(1)=[];
fig_oil=figure;
subplot(3,1,1); plot(Y(id_oil)+M(id_oil)/12+D(id_oil)/365,([oil_vec(id_oil,1),model_oil_vec(id_oil,1)])); title('1m oil futures');
subplot(3,1,2); plot(Y(id_oil)+M(id_oil)/12+D(id_oil)/365,([oil_vec(id_oil,2),model_oil_vec(id_oil,2)])); title('3m oil futures');
subplot(3,1,3); plot(Y(id_oil)+M(id_oil)/12+D(id_oil)/365,([oil_vec(id_oil,3),model_oil_vec(id_oil,3)])); title('12m oil futures');
% print(fig_oil,'-depsc2','..\results\oil_FRBAoil.eps');


fig_IE=figure;
plot(Y+M/12+D/365,IE,'LineWidth',3); title('Inflation Expectations','FontSize',12); 
legend('1-year','10-year');

% print(fig_IE,'-depsc2','..\results\IE_FRBAoil.eps');

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

% print(fig_allIE,'-depsc2','..\results\allIE_FRBAoil.eps');


% % %%%%%%%%%%%%%%%%%%%%%%%%
% % % Forward Term Premium %
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % load ('..\data\allyield.txt');
% % allyield(allyield(:,1)<1990,:)=[];
% % 
% % yields=allyield(:,[10,15])/100;
% % 
% % mats=[5 10]';
% % T=size(yields,1);
% % prices=-(ones(T,1)*mats').*yields;
% % forwards = prices(:,1)-prices(:,2);
% % 
% % 
% % % Compute the nominal yield factor loadings
% % [ay_tran,by_tran,cy_tran] = YieldFacLoad_ODE3(5,tmp_k,tmp_K,tmp_h,tmp_H,rho0,[rho1;rhod;rhos],rhov,KAPPAv_rn,SIGMAv^2);
% % ay=ay_tran; by=by_tran(1:Nfac,:); cy=cy_tran; dy=by_tran(Nfac+1,:); ey=by_tran(Nfac+2,:);
% % 
% % ay=-ay./5; ay=ay';
% % by=-by./repmat(5,Nfac,1); by=by';
% % cy=-cy./5; cy=cy';
% % dy=-dy./5; dy=dy';
% % ey=-ey./5; ey=ey';
% % 
% % FTP=zeros(T,1);
% % for i=1:T
% %     ForecastHor=5;
% %     xt=state(i,2:4)'; vt=state(i,5); dt=state(i,6); st=state(i,7);
% % 
% %    af=ay+by*(eye(Nfac)-expm(-KAPPA*ForecastHor))*theta ...
% %          + cy*(1-exp(-KAPPAv*ForecastHor))*thetav ...
% %          + dy*(1-exp(-KAPPAd*ForecastHor))*thetad ...
% %          + ey*(1/phis*(exp(phis*ForecastHor)-1)*(phi0+phi1'*theta+phiv*thetav+phid*thetad)) ...
% %          - ey*phi1'*inv(KAPPA+phis*eye(Nfac))*(exp(phis*ForecastHor)*eye(Nfac)-expm(-KAPPA*ForecastHor))*theta ...
% %          - ey*phiv*inv(KAPPAv+phis)*(exp(phis*ForecastHor)-exp(-KAPPAv*ForecastHor))*thetav ...
% %          - ey*phid*inv(KAPPAd+phis)*(exp(phis*ForecastHor)-exp(-KAPPAd*ForecastHor))*thetad;         
% %   bf=by*expm(-KAPPA*ForecastHor) ...
% %          +ey*phi1'*inv(KAPPA+phis*eye(Nfac))*(exp(phis*ForecastHor)*eye(Nfac)-expm(-KAPPA*ForecastHor));
% %   cf=cy*exp(-KAPPAv*ForecastHor) ...
% %          +ey*phiv*inv(KAPPAv+phiv)*(exp(phis*ForecastHor)-exp(-KAPPAv*ForecastHor));
% %   df=dy*exp(-KAPPAd*ForecastHor) ...
% %          +ey*phid*inv(KAPPAd+phid)*(exp(phis*ForecastHor)-exp(-KAPPAd*ForecastHor));
% %   ef=ey*exp(phis*ForecastHor);
% % 
% %     
% %     tmp=af+bf*xt+cf*vt+df*dt+ef*st;
% %     FTP(i)=forwards(i)/5-tmp;      
% % end
% % plot([FTP,state(:,5)])
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP,[ones(T,1),state(:,5)],12,0)
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP,[ones(T,1),state(:,2:4)],12,0)
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP(id_floors),[ones(length(id_floors),1),state(id_floors,5)],12,0)
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP(id_floors),[ones(length(id_floors),1),state(id_floors,2:4)],12,0)
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % % Forward Term Premium (REAL) %
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % load ('..\data\alltips.txt');
% % alltips(alltips(:,1)<1990,:)=[];
% % 
% % yields=alltips(:,[7,12])/100;
% % 
% % mats=[5 10]';
% % T=size(yields,1);
% % prices=-(ones(T,1)*mats').*yields;
% % forwards = prices(:,1)-prices(:,2);
% % 
% % 
% % % Compute the nominal yield factor loadings
% % [ay_R,by_R,cy_R,dy_R] = YieldFacLoad_ODE4(5,KAPPAtheta_rn_R,KAPPA_rn_R,KAPPAv*thetav,SIGMA*SIGMA', ...
% %                                           rho0_R,rho1_R,rhov_R,rhod_R,KAPPAv_rn,SIGMAv^2,KAPPAd_rn*thetad_rn,KAPPAd_rn,SIGMAd^2);
% % 
% % ay_R=-ay_R./5; ay_R=ay_R';
% % by_R=-by_R./5; by_R=by_R';
% % cy_R=-cy_R./5; cy_R=cy_R';
% % dy_R=-dy_R./5; dy_R=dy_R';
% % 
% % 
% % FTP_R=zeros(T,1);
% % for i=1:T
% %     xt=state(i,2:4)'; vt=state(i,5); dt=state(i,6);
% %     tmp=ay_R+by_R*(theta+expm(-KAPPA*5)*(xt-theta)) ...
% %             +cy_R*(thetav+exp(-KAPPAv*5)*(vt-thetav)) ...
% %             +dy_R*(thetad+exp(-KAPPAd*5)*(dt-thetad));
% %     FTP_R(i)=forwards(i)/5-tmp;
% % end
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP_R,[ones(T,1),state(:,5)],12,0)
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP_R,[ones(T,1),state(:,2:4)],12,0)
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP_R(id_floors),[ones(length(id_floors),1),state(id_floors,5)],12,0)
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(FTP_R(id_floors),[ones(length(id_floors),1),state(id_floors,2:4)],12,0)
% % 
% % % figure; plot(Y+((M-1)*30+D)/365,[FTP,FTP_R]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bond Return Predictability %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('..\data\allyield.txt');
allyield(allyield(:,1)<1990,:)=[];
allyield(allyield(:,1)>2014,:)=[];

yields=allyield(:,6:10)/100; 
mats=[1 2 3 4 5]';
T=size(yields,1);

prices=-(ones(T,1)*mats').*yields;
forwards = prices(:,1:4)-prices(:,2:5);
fs = forwards-yields(:,1)*ones(1,4);

% hprx(t) is the holding period return over last year

hpr = prices(53:T,1:4)-prices(1:T-52,2:5);
hprx = hpr - yields(1:T-52,1)*ones(1,4);
hpr = [zeros(52,1)*ones(1,4); hpr];     % pads out the initial values with zeros so same length as other series       
hprx = [zeros(52,1)*ones(1,4); hprx];

% capitalized variables do not follow the convention of being padded out with initial zeros
% instead, HPRX starts 12 months later, so the first HPRX is in 65 while the first FS, FT, YT is 1964. 
% These are set up so you can regress HPRX, AHPRX on YT, FT, etc. directly
% with no subscripting. They also include a column of ones for use in
% regressions. 

HPRX = 100*hprx(53:T,:);
AHPRX = mean(HPRX')';

% CP factor %
FT = [ones(T-52,1) yields(1:T-52,1) forwards(1:T-52,:)]; % yeilds and forwards
[gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,FT,12,0); % std errors using HH
disp(' gammas, ols se');
disp([ gammas'; olsse'])
disp(' R2 ');
disp(R2hump);

% % CP=FT_all*gammas;
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(CP,[ones(T,1),state(:,6)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(CP,[ones(T,1),state(:,5)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(CP,[ones(T,1),state(:,4)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(CP,[ones(T,1),state(:,3)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(CP,[ones(T,1),state(:,2)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),CP(1:end-52)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump


% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,5)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,6)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,5:6)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,2)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,3)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,4)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),state(1:end-52,2:6)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),CP(1:end-52),state(1:end-52,2:6)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump
% % 
% % i=2
% % [gammas,olsse,R2hump,R2humpadj,v] = olsgmm(AHPRX,[ones(T-52,1),CP(1:end-52),state(1:end-52,i)],12,0); % std errors using HH
% % gammas'
% % (gammas./olsse)'
% % R2hump


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % In-Sample Inflation Fitting %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % hor_is=[1 6 12];    % forecasting horizon
% % Nhor_is=length(hor_is);
% % 
% % Np=length(id_p);
% % phat_is=NaN(Np,Nhor_is);
% % phat_is_RW=NaN(Np,Nhor_is);
% % 
% % p_is_rmse=zeros(Nhor_is,1);
% % p_is_RW_rmse=zeros(Nhor_is,1);
% % 
% % for j=1:Nhor_is
% %     p_prev=exp(p_vec(id_p(1)));
% %     hor=hor_is(j);
% %     
% %     for i=(hor+1):Np
% %         MAT=(id_p(i)-id_p(i-hor))*Dt;
% %                 
% %         ax=(eye(Nfac)-inv(KAPPA*MAT)*(eye(Nfac)-expm(-KAPPA*MAT)))*theta;
% %         bx=(eye(Nfac)-expm(-KAPPA'*MAT))*inv(KAPPA'*MAT);
% %         ad=(1-inv(KAPPAd*MAT)*(1-exp(-KAPPAd*MAT)))*thetad;
% %         bd=(1-exp(-KAPPAd'*MAT))*inv(KAPPAd'*MAT);
% % 
% %         W=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*MAT))/MAT;
% %         Wd=inv(KAPPAd)*(1-exp(-KAPPAd*MAT))/MAT;
% %         Ws=inv(phis)*(exp(phis*MAT)-1)/MAT;
% % 
% %         bsx=phi1'*inv(KAPPA+phis*eye(Nfac))*(Ws*eye(Nfac)-W);
% %         bsd=phid*inv(KAPPAd+phis)*(Ws-Wd);        
% %         bss=Ws;
% %         as=1/phis*(Ws-1)*(phi0+phi1'*theta+phid*thetad) ...
% %              -bsx*theta-bsd*thetad; 
% % 
% %         aI=rho0+rho1'*ax+rhod*ad+rhos*as;
% %         bI=rho1'*bx+rhos*bsx; bI=bI';
% %         cI2=rhod*bd+rhos*bsd;
% %         cI3=rhos*bss;
% % 
% %         phat_is(i,j)=p_prev*exp((aI+state(id_p(i-hor),2:4)*bI+state(id_p(i-hor),5)*cI2+state(id_p(i-hor),6)*cI3)*MAT);
% % 
% %         phat_is_RW(i,j)=p_prev;
% % 
% %         p_prev=exp(p_vec(id_p(i-hor+1)));
% %     end
% %     fe_is=(phat_is((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
% %     fe_is_RW=(phat_is_RW((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
% %     
% %     p_is_rmse(j)=sqrt(mean((fe_is.^2)));
% %     p_is_RW_rmse(j)=sqrt(mean(fe_is_RW.^2));
% % 
% % end
% % [p_is_rmse,p_is_RW_rmse]
% % 
% % figure;
% % for j=1:3
% %     subplot(3,1,j);
% %     plot([exp(p_vec(id_p((1+hor_is(j)):end))),phat_is((1+hor_is(j)):end,j),phat_is_RW((1+hor_is(j)):end,j)])
% %     legend('data','model','rw');
% % end
% % 
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Out-of-sample (OOS) Inflation Fitting %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % id_lastobs=find(Y==2013 & M==12 & D==18);
% % state_lastobs=state(id_lastobs,:);
% % p_lastobs=exp(p_vec(id_lastobs));
% % 
% % startdaten=datenum(2013,12,18);
% % enddaten=datenum(2015,12,31);
% % 
% % data_FRBA; 
% % [Y,M,D]=datevec(mydate);
% % 
% % Dt=7/365;
% % % hor_oos=[1 6 12];    % forecasting horizon in months
% % % Nhor_oos=length(hor_oos);
% % 
% % Np=length(id_p);
% % phat_oos=p_lastobs*ones(Np,1);
% % phat_oos_RW=p_lastobs*ones(Np,1);
% % 
% % 
% % for i=2:Np
% %     MAT=(id_p(i)-id_p(1))*Dt;
% %     
% %     ax=(eye(Nfac)-inv(KAPPA*MAT)*(eye(Nfac)-expm(-KAPPA*MAT)))*theta;
% %     bx=(eye(Nfac)-expm(-KAPPA'*MAT))*inv(KAPPA'*MAT);
% %     ad=(1-inv(KAPPAd*MAT)*(1-exp(-KAPPAd*MAT)))*thetad;
% %     bd=(1-exp(-KAPPAd'*MAT))*inv(KAPPAd'*MAT);
% %     
% %     W=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*MAT))/MAT;
% %     Wd=inv(KAPPAd)*(1-exp(-KAPPAd*MAT))/MAT;
% %     Ws=inv(phis)*(exp(phis*MAT)-1)/MAT;
% % 
% %     bsx=phi1'*inv(KAPPA+phis*eye(Nfac))*(Ws*eye(Nfac)-W);
% %     bsd=phid*inv(KAPPAd+phis)*(Ws-Wd);        
% %     bss=Ws;
% %     as=1/phis*(Ws-1)*(phi0+phi1'*theta+phid*thetad) ...
% %          -bsx*theta-bsd*thetad; 
% %         
% %     aI=rho0+rho1'*ax+rhod*ad+rhos*as;
% %     bI=rho1'*bx+rhos*bsx; bI=bI';
% %     cI2=rhod*bd+rhos*bsd;
% %     cI3=rhos*bss;
% % 
% %     phat_oos(i)=p_prev*exp((aI+state_lastobs(:,2:4)*bI+state_lastobs(:,5)*cI2+state_lastobs(:,6)*cI3)*MAT);
% % 
% % end
% % fe_oos=(phat_oos-exp(p_vec(id_p)))./exp(p_vec(id_p));
% % fe_oos_RW=(phat_oos_RW-exp(p_vec(id_p)))./exp(p_vec(id_p));
% % 
% % p_oos_rmse=sqrt(mean((fe_oos.^2)));
% % p_oos_RW_rmse=sqrt(mean(fe_oos_RW.^2));
% % 
% % [p_oos_rmse,p_oos_RW_rmse]
% % 
% % figure;
% % plot([exp(p_vec(id_p)),phat_oos,phat_oos_RW])
% % legend('data','model','rw');
% % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample Inflation Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax=(eye(Nfac)-inv(KAPPA*Dt)*(eye(Nfac)-expm(-KAPPA*Dt)))*theta;
bx=(eye(Nfac)-expm(-KAPPA'*Dt))*inv(KAPPA'*Dt);
ad=(1-inv(KAPPAd*Dt)*(1-exp(-KAPPAd*Dt)))*thetad;
bd=(1-exp(-KAPPAd'*Dt))*inv(KAPPAd'*Dt);

W=inv(KAPPA)*(eye(Nfac)-expm(-KAPPA*Dt))/Dt;
Wd=inv(KAPPAd)*(1-exp(-KAPPAd*Dt))/Dt;
Ws=inv(phis)*(exp(phis*Dt)-1)/Dt;

bsx=phi1'*inv(KAPPA+phis*eye(Nfac))*(Ws*eye(Nfac)-W);
bsd=phid*inv(KAPPAd+phis)*(Ws-Wd);        
bss=Ws;
as=1/phis*(Ws-1)*(phi0+phi1'*theta+phid*thetad) ...
     -bsx*theta-bsd*thetad; 

aI=rho0_pi+rho1_pi'*ax+rhod_pi*ad+rhos_pi*as;
bI=rho1_pi'*bx+rhos_pi*bsx; bI=bI';
cI2=rhod_pi*bd+rhos_pi*bsd;
cI3=rhos_pi*bss;
 
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

        phat_is(i,j)=p_prev*exp((aI*length(id_forecast)+sum(state(id_forecast,2:4))*bI ...
                    +sum(state(id_forecast,5))*cI2+sum(state(id_forecast,6))*cI3)*Dt);

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
% Out-of-sample (OOS) Inflation Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startdaten=datenum(1990,1,1);
enddaten=datenum(2015,12,31);
data_FRBA; 
[Y,M,D]=datevec(mydate);
[logL,logL_vec,state,model_oil_vec]= logL_FRBAoil(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,oil_vec,MATgrid_oil,notmissing_oil);

id_lastobs=find(Y==2013 & M==12 & D==18);
p_vec(1:(id_lastobs-1))=[];
state(1:(id_lastobs-1),:)=[];


hor_oos=[1 6];    % forecasting horizon in months
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
        phat_oos(i)=p_prev*exp((aI*length(id_forecast)+sum(state(id_forecast,2:4))*bI ...
                    +sum(state(id_forecast,5))*cI2+sum(state(id_forecast,6))*cI3)*Dt);

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



