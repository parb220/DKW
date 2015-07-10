% function [logL,state,IE,IRP,FTP,FTP_R]=output_DKWl_os_9013

% close all; clear; clc;

% startdaten=datenum(1990,1,1);
% enddaten=datenum(2013,12,31);
% data_FRBA; 
% [Y,M,D]=datevec(mydate);
% load '.\results\pvar_DKW_9013.mat';
% id_lastobs=find(Y==2013 & M==12 & D==18);

is_year_beg=1990; % beginning year of in-sample period
is_year_end=2012; % ending year of in-sample period

startdaten=datenum(is_year_beg,1,1);
enddaten=datenum(is_year_end,12,31);
load '.\results\pvar_DKWl_os_9012.mat';
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
sigq=NaN(Nfac,1);
sigqx=NaN;

delta_p=NaN;
delta_tips=NaN(Nytips,1);
delta_options=NaN(length(MATgrid_options),1);
delta_swap=NaN(Nytips,1);   % same number of swaps as TIPS

KAPPA_L=NaN;       
SIGMA_L=NaN;
theta_L=NaN;   
rho1_L=NaN(Nfac,1);
rhoL_L=NaN;
lambda0_L=NaN;
SIGMAlambda1_L=NaN;

% parameter normalization
KAPPA=diag(NaN(Nfac,1));
SIGMA(1,1)=0.01;
SIGMA(1,2)=0;
SIGMA(1,3)=0;
SIGMA(2,2)=0.01;
SIGMA(2,3)=0;
SIGMA(3,3)=0.01;

theta=zeros(Nfac,1);

SIGMA_L=0.01;

delta_bcfLT=0.0075;

delta_p=0;

paras_vec=[reshape(KAPPA,Nfac^2,1);reshape(SIGMA,Nfac^2,1);theta;
      rho0;rho1;lambda0;reshape(SIGMAlambda1,Nfac^2,1);
      delta_y;delta_bcf;delta_bcfLT;rho0_pi;rho1_pi;sigq;
      sigqx;delta_p;delta_tips;delta_options;delta_swap;KAPPA_L;SIGMA_L;theta_L;
      rho1_L;rhoL_L;lambda0_L;SIGMAlambda1_L];


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

delta_swap=pall(3*Nfac^2+5*Nfac+8+Ny+Nytips+length(MATgrid_options):3*Nfac^2+5*Nfac+7+Ny+2*Nytips+length(MATgrid_options));
KAPPA_L=pall(3*Nfac^2+5*Nfac+8+Ny+2*Nytips+length(MATgrid_options));
SIGMA_L=pall(3*Nfac^2+5*Nfac+9+Ny+2*Nytips+length(MATgrid_options));
theta_L=pall(3*Nfac^2+5*Nfac+10+Ny+2*Nytips+length(MATgrid_options));
rho1_L=pall(3*Nfac^2+5*Nfac+11+Ny+2*Nytips+length(MATgrid_options):3*Nfac^2+6*Nfac+10+Ny+2*Nytips+length(MATgrid_options));
rhoL_L=pall(3*Nfac^2+6*Nfac+11+Ny+2*Nytips+length(MATgrid_options));
lambda0_L=pall(3*Nfac^2+6*Nfac+12+Ny+2*Nytips+length(MATgrid_options));
SIGMAlambda1_L=pall(3*Nfac^2+6*Nfac+13+Ny+2*Nytips+length(MATgrid_options));


lambda1=inv(SIGMA)*SIGMAlambda1;
lambda1_L=inv(SIGMA_L)*SIGMAlambda1_L;

lambda1=inv(SIGMA)*SIGMAlambda1;
lambda1_L=inv(SIGMA_L)*SIGMAlambda1_L;

[logL,logL_vec,state]= logL_DKWl(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT);


Np=length(pvar);
% % rowLabels = {'row 1','row 2','row 3','row 4','row 5','row 6','row 7','row 8','row 9','row 10', ...
% %              'row 11','row 12','row 13','row 14','row 15','row 16','row 17','row 18','row 19','row 20', ...
% %              'row 21','row 22','row 23','row 24','row 25','row 26','row 27','row 28','row 29','row 30', ...
% %              'row 31','row 32','row 33','row 34','row 35','row 36','row 37','row 38','row 39','row 40', ...
% %              'row 41','row 42','row 43','row 44','row 45','row 46','row 47','row 48','row 49','row 50','row 51'};
% % columnLabels = {'pvar'};
% % 
% % pvar_mat=NaN(51,1);
% % pvar_mat(1:3,:)=[pvar(1:3)]; % [KAPPA(1,1);KAPPA(2,2);KAPPA(3,3)];
% % pvar_mat(4:6,:)=100*[pvar(4:6)]; % [SIGMA(2,1);SIGMA(3,1);SIGMA(3,2)];
% % pvar_mat(7:10,:)=[pvar(7:10)]; % [rho0;rho1];
% % pvar_mat(11,:)=[NaN]; % rhod
% % pvar_mat(12:23,:)=[pvar(11:22)]; % [lambda0;SIGMAlambda1(:)];
% % pvar_mat(24:27,:)=[pvar(32:35)]; % [rho0_pi; rho1_pi];
% % pvar_mat(28,:)=[NaN]; % rhod_pi
% % pvar_mat(29:32,:)=100*[pvar(36:39)]; % [sigq;sigqx];
% % pvar_mat(33:38,:)=NaN(6,1); % [KAPPAv;thetav;SIGMAv;gamv;gamvx;rho];
% % pvar_mat(39:45,:)=100*[pvar(23:29)]; % delta_y
% % pvar_mat(46:48,:)=100*[pvar(40:42)]; % delta_tips;
% % pvar_mat(49:50,:)=100*[pvar(30:31)]; % [delta_bcf];
% % pvar_mat(51,:)=[NaN]; %[delta_caps];
% % 
% % % matrix2latex(pvar_mat, '.\results\pvar_DKW_9013.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');



% % % parameters: real pricing kernel (derived from above)
rho0_R=rho0-rho0_pi-1/2*(sigq'*sigq+sigqx^2)+lambda0'*sigq;
rho1_R=rho1-rho1_pi+lambda1'*sigq;
lambda0_R=lambda0-sigq;
SIGMAlambda1_R=SIGMAlambda1;

KAPPA_R_rn=KAPPA + SIGMAlambda1_R;
KAPPAtheta_R_rn=KAPPA*theta - SIGMA*lambda0_R;

% Compute factor loadings for LIQUIDITY component of TIPS yields
KAPPA_L_rn=KAPPA_L + SIGMAlambda1_L;
KAPPAtheta_L_rn=KAPPA_L*theta_L - SIGMA_L*lambda0_L;


%%%%%%%%%%%%%%%%%%%%%%%
% Yield Decomposition %
%%%%%%%%%%%%%%%%%%%%%%%

horLIST=[1 10];
IE=zeros(T,length(horLIST));
IRP=zeros(T,length(horLIST));
LRP=zeros(T,length(horLIST));

for i=1:length(horLIST);
    MAT=horLIST(i);

    [ay1,by1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MAT);    
    
    [Ay_R,By_R]=YieldFacLoad_ODE(MAT,KAPPAtheta_R_rn,KAPPA_R_rn,SIGMA*SIGMA',rho0_R,rho1_R);
    ay_R=-Ay_R./MAT; ay_R=ay_R'; 
    by_R=-By_R./repmat(MAT,Nfac,1); by_R=by_R';
    ay_R1=ay_R; by_R1=by_R;


    % Compute factor loadings for REAL component of TIPS yields
    [Ay_TR,By_TR]=YieldFacLoad_ODE(MAT,KAPPAtheta_R_rn,KAPPA_R_rn,SIGMA*SIGMA',rho0_R,rho1_R+rho1_L);
    ay_TR=-Ay_TR./MAT; ay_TR=ay_TR';
    by_TR=-By_TR./repmat(MAT,Nfac,1); by_TR=by_TR';
    ay_TR1=ay_TR; by_TR1=by_TR;

    [Ay_L,By_L]=YieldFacLoad_ODE(MAT,KAPPAtheta_L_rn,KAPPA_L_rn,SIGMA_L*SIGMA_L',0,rhoL_L);
    ay_L=-Ay_L./MAT; ay_L=ay_L';
    by_L=-By_L./MAT; by_L=by_L';
    ay_L1=ay_L; by_L1=by_L;

    aI1=rho0_pi+rho1_pi'*(eye(Nfac)-1/MAT*inv(-KAPPA)*(expm(-KAPPA*MAT)-eye(Nfac)))*theta;
    bI1=1/MAT*inv(-KAPPA')*(expm(-KAPPA'*MAT)-eye(Nfac))*rho1_pi; 
    bI1=bI1';

    IE(:,i)=aI1+state(:,2:4)*bI1';
    IRP(:,i)=(ay1-ay_R1-aI1)+state(:,2:4)*(by1-by_R1-bI1)';
    LRP(:,i)=(ay_TR1+state(:,2:4)*by_TR1')+(ay_L1+state(:,5)*by_L1')-(ay_R1+state(:,2:4)*by_R1');
    
    if MAT==10
        BEI10_act_tips=y_vec(:,end)-ytips_vec(:,end);
        BEI10_model_tips=(ay1+state(:,2:4)*by1')-((ay_TR1+state(:,2:4)*by_TR1')+(ay_L1+state(:,5)*by_L1'));
        BEI10_model_true=(ay1+state(:,2:4)*by1')-(ay_R1+state(:,2:4)*by_R1');
    end
end

fig_allBEI=figure;
plot(Y(id_tips)+M(id_tips)/12+D(id_tips)/365,BEI10_act_tips(id_tips),'r'); hold on;
plot(Y(id_tips)+M(id_tips)/12+D(id_tips)/365,BEI10_model_tips(id_tips),'k'); hold on;
plot(Y+M/12+D/365,BEI10_model_true,'b'); hold off;
legend('Actual TIPS','Model TIPS','Model True');


fig_state=figure; 
subplot(2,2,1); plot(Y+M/12+D/365,state(:,2),'LineWidth',4); title('DKW: state variable x1','FontSize',12);
subplot(2,2,2); plot(Y+M/12+D/365,state(:,3),'LineWidth',4); title('DKW: state variable x2','FontSize',12);
subplot(2,2,3); plot(Y+M/12+D/365,state(:,4),'LineWidth',4); title('DKW: state variable x3','FontSize',12); ylim([-0.2 0.2]);

fig_IE=figure;
plot(Y+M/12+D/365,IE,'LineWidth',3); title('Inflation Expectations','FontSize',12); 
legend('1-year','10-year');

fig_IRP=figure;
plot(Y+M/12+D/365,IRP,'LineWidth',3); title('Inflation Risk Premiums','FontSize',12); 
legend('1-year','10-year');

fig_LRP=figure;
plot(Y+M/12+D/365,LRP,'LineWidth',3); title('Liquidity Risk Premiums','FontSize',12); 
legend('1-year','10-year');



% % %%%%%%%%%%%%%%%%%%%%%%%%
% % % Forward Term Premium %
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % load ('..\data\allyield.txt');
% % allyield(allyield(:,1)<1990,:)=[];
% % 
% % yields=allyield(:,[8,13])/100;
% % 
% % mats=[5 10]';
% % T=size(yields,1);
% % prices=-(ones(T,1)*mats').*yields;
% % forwards = prices(:,1)-prices(:,2);
% % 
% % 
% % % Compute the nominal yield factor loadings
% % [ay_DKW,by_DKW]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,5);
% % 
% % FTP=zeros(T,1);
% % for i=1:T  
% %     xt_DKW=state(i,2:4)';
% %     tmp_DKW=ay_DKW+by_DKW*(theta+expm(-KAPPA*5)*(xt_DKW-theta));
% %     FTP(i)=forwards(i)/5-tmp_DKW;
% %     
% % end
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
% % rho0_R=rho0-rho0_pi-1/2*(sigq'*sigq+sigqx^2)+lambda0'*sigq;
% % rho1_R=rho1-rho1_pi+lambda1'*sigq;
% % lambda0_R=lambda0-sigq;
% % SIGMAlambda1_R=SIGMAlambda1;
% % 
% % [ay_R,by_R]=YieldFacLoad(KAPPA,SIGMA,theta,rho0_R,rho1_R,lambda0_R,SIGMAlambda1_R,5);
% % 
% % 
% % FTP_R=zeros(T,1);
% % for i=1:T
% %     xt=state(i,2:4)'; 
% %     tmp=ay_R+by_R*(theta+expm(-KAPPA*5)*(xt-theta));
% %     FTP_R(i)=forwards(i)/5-tmp;
% % end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % In-Sample Inflation Fitting %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dt=7/365;
% hor_is=[1 6 12];    % forecasting horizon in months
% Nhor_is=length(hor_is);
% 
% Np=length(id_p);
% phat_is=NaN(Np,Nhor_is);
% phat_is_RW=NaN(Np,Nhor_is);
% 
% p_is_rmse=zeros(Nhor_is,1);
% p_is_RW_rmse=zeros(Nhor_is,1);
% 
% for j=1:Nhor_is
%     p_prev=exp(p_vec(id_p(1)));
% 
%     hor=hor_is(j);
%     
%     for i=(hor+1):Np
%         MAT=(id_p(i)-id_p(i-hor))*Dt;
%         ax=(eye(Nfac)-inv(KAPPA*MAT)*(eye(Nfac)-expm(-KAPPA*MAT)))*theta;
%         bx=(eye(Nfac)-expm(-KAPPA'*MAT))*inv(KAPPA'*MAT);
% 
%         aI1=rho0_pi+rho1_pi'*ax;
%         bI1=bx'*rho1_pi; bI1=bI1';
%         
%         phat_is(i,j)=p_prev*exp((aI1+state(id_p(i-hor),2:4)*bI1')*MAT);
% 
%         phat_is_RW(i,j)=p_prev;
% 
%         p_prev=exp(p_vec(id_p(i-hor+1)));
%     end
%     fe_is=(phat_is((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
%     fe_is_RW=(phat_is_RW((hor+1):Np,j)-exp(p_vec(id_p((hor+1):Np))))./exp(p_vec(id_p((hor+1):Np)));
%     
%     p_is_rmse(j)=sqrt(mean((fe_is.^2)));
%     p_is_RW_rmse(j)=sqrt(mean(fe_is_RW.^2));
% 
% end
% [p_is_rmse,p_is_RW_rmse]
% 
% 
% figure;
% for j=1:3
%     subplot(3,1,j);
%     plot([exp(p_vec(id_p((1+hor_is(j)):end))),phat_is((1+hor_is(j)):end,j),phat_is_RW((1+hor_is(j)):end,j)])
%     legend('data','model','rw');
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Out-of-sample (OOS) Inflation Fitting %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% id_lastobs=find(Y==2013 & M==12 & D==18);
% state_lastobs=state(id_lastobs,:);
% p_lastobs=exp(p_vec(id_lastobs));
% 
% startdaten=datenum(2013,12,18);
% enddaten=datenum(2015,12,31);
% 
% data_FRBA; 
% [Y,M,D]=datevec(mydate);
% 
% Dt=7/365;
% 
% Np=length(id_p);
% phat_oos=p_lastobs*ones(Np,1);
% phat_oos_RW=p_lastobs*ones(Np,1);
% 
% 
% for i=2:Np
%     MAT=(id_p(i)-id_p(1))*Dt;
%     ax=(eye(Nfac)-inv(KAPPA*MAT)*(eye(Nfac)-expm(-KAPPA*MAT)))*theta;
%     bx=(eye(Nfac)-expm(-KAPPA'*MAT))*inv(KAPPA'*MAT);
% 
%     aI1=rho0_pi+rho1_pi'*ax;
%     bI1=bx'*rho1_pi; bI1=bI1';
% 
%     phat_oos(i)=p_prev*exp((aI1+state_lastobs(:,2:4)*bI1')*MAT);
% 
% end
% fe_oos=(phat_oos-exp(p_vec(id_p)))./exp(p_vec(id_p));
% fe_oos_RW=(phat_oos_RW-exp(p_vec(id_p)))./exp(p_vec(id_p));
% 
% p_oos_rmse=sqrt(mean((fe_oos.^2)));
% p_oos_RW_rmse=sqrt(mean(fe_oos_RW.^2));
% 
% [p_oos_rmse,p_oos_RW_rmse]
% 
% figure;
% plot([exp(p_vec(id_p)),phat_oos,phat_oos_RW])
% legend('data','model','rw');

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
        
    [ay1,by1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MAT);    
    
    for j=1:Nhor_is
        hor=hor_is(j);

        yNhat_is=NaN(T,1);
        yNhat_is_RW=NaN(T,1);

        for i=(hor+1):T
            yNhat_is(i)=ay1+by1*(theta+expm(-KAPPA*hor*Dt)*(state(i-hor,2:4)'-theta));
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
        
    [ay_R1,by_R1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0_R,rho1_R,lambda0_R,SIGMAlambda1_R,MAT);
    
    for j=1:Nhor_is
        hor=hor_is(j);

        yRhat_is=NaN(T,1);
        yRhat_is_RW=NaN(T,1);

        for i=(hor+id_alltips_start):T
            yRhat_is(i)=ay_R1+by_R1*(theta+expm(-KAPPA*hor*Dt)*(state(i-hor,2:4)'-theta));
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

aI1=rho0_pi+rho1_pi'*ax;
bI1=bx'*rho1_pi; bI1=bI1';


hor_is=[1 6 12];    % forecasting horizon in months
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
        
        phat_is(i,j)=p_prev*exp((aI1*length(id_forecast)+sum(state(id_forecast,2:4))*bI1')*Dt);
        
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
[logL,logL_vec,state]= logL_DKW(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT);


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
        
    [ay1,by1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0,rho1,lambda0,SIGMAlambda1,MAT);
        
    for j=1:Nhor_oos
        hor=hor_oos(j);
    
        yNhat_oos=NaN(T_oos,1);
        yNhat_oos_RW=NaN(T_oos,1);

        for i=(hor+1):T_oos
            yNhat_oos(i)=ay1+by1*(theta+expm(-KAPPA*hor*Dt)*(state_oos(i-hor,2:4)'-theta));
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
        
    [ay_R1,by_R1]=YieldFacLoad(KAPPA,SIGMA,theta,rho0_R,rho1_R,lambda0_R,SIGMAlambda1_R,MAT);
        
    for j=1:Nhor_oos
        hor=hor_oos(j);
    
        yRhat_oos=NaN(T_oos,1);
        yRhat_oos_RW=NaN(T_oos,1);

        for i=(hor+1):T_oos
            yRhat_oos(i)=ay_R1+by_R1*(theta+expm(-KAPPA*hor*Dt)*(state_oos(i-hor,2:4)'-theta));
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
        id_forecast=id_p_oos(i-hor):(id_p_oos(i)-1);
        phat_oos(i)=p_prev*exp((aI1*length(id_forecast)+sum(state(id_forecast,2:4))*bI1')*Dt);
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


