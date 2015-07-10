
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
load '.\results\pvar_DKW_nominal_9012.mat';
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

rho0_pi=NaN;
rho1_pi=NaN(Nfac,1);
sigq=NaN(Nfac,1);
sigqx=NaN;

delta_p=NaN;

% parameter normalization
KAPPA=diag(NaN(Nfac,1));
SIGMA(1,1)=0.01;
SIGMA(1,2)=0;
SIGMA(1,3)=0;
SIGMA(2,2)=0.01;
SIGMA(2,3)=0;
SIGMA(3,3)=0.01;

theta=zeros(Nfac,1);
delta_p=0;
    
paras_vec=[reshape(KAPPA,Nfac^2,1);reshape(SIGMA,Nfac^2,1);theta;
      rho0;rho1;lambda0;reshape(SIGMAlambda1,Nfac^2,1);delta_y;
      rho0_pi;rho1_pi;sigq;sigqx;delta_p];

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

[logL,logL_vec,state]= logL_DKW_nominal(pvar,paras_vec,pidx,Nfac,Dt,y_vec,MATgrid,p_vec,notmissing_p);


Np=length(pvar);
rowLabels = {'row 1','row 2','row 3','row 4','row 5','row 6','row 7','row 8','row 9','row 10', ...
             'row 11','row 12','row 13','row 14','row 15','row 16','row 17','row 18','row 19','row 20', ...
             'row 21','row 22','row 23','row 24','row 25','row 26','row 27','row 28','row 29','row 30', ...
             'row 31','row 32','row 33','row 34','row 35','row 36','row 37','row 38','row 39','row 40', ...
             'row 41','row 42','row 43','row 44','row 45','row 46','row 47','row 48','row 49','row 50'};
columnLabels = {'pvar'};

pvar_mat=NaN(50,1);
pvar_mat(1:3,:)=[KAPPA(1,1);KAPPA(2,2);KAPPA(3,3)]; % [KAPPA(1,1);KAPPA(2,2);KAPPA(3,3)];
pvar_mat(4:6,:)=100*[SIGMA(2,1);SIGMA(3,1);SIGMA(3,2)]; % [SIGMA(2,1);SIGMA(3,1);SIGMA(3,2)];
pvar_mat(7:10,:)=[rho0;rho1]; % [rho0;rho1];
pvar_mat(11,:)=[NaN]; % rhov
pvar_mat(12:23,:)=[lambda0;SIGMAlambda1(:)]; % [lambda0;SIGMAlambda1(:)];
pvar_mat(24:27,:)=[rho0_pi; rho1_pi]; % [rho0_pi; rho1_pi];
pvar_mat(28,:)=[NaN]; % rhov_pi
pvar_mat(29:32,:)=100*[sigq;sigqx]; % [sigq;sigqx];
pvar_mat(33:38,:)=NaN(6,1); % [KAPPAv;thetav;SIGMAv;gamv;gamvx;rho];
pvar_mat(39:45,:)=100*[delta_y]; % delta_y
pvar_mat(46:48,:)=100*[NaN;NaN;NaN]; % delta_tips;
pvar_mat(49:50,:)=100*[NaN;NaN]; % [delta_bcf];

matrix2latex(pvar_mat, '.\results\pvar_DKW_nominal_9012.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-Sample Nominal Yields Fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT_allyields=[0.25 0.5 1:30];
MAT_is=[0.5 1:10]; % bond maturity in years (6m/2y/10y)
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
[MAT_list,HOR_list,yN_is_RW_rmse(:)*100,yN_is_rmse(:)*100]

rowLabels = {'row 1','row 2','row 3','row 4','row 5','row 6','row 7','row 8','row 9','row 10', ...
             'row 11','row 12','row 13','row 14','row 15','row 16','row 17','row 18','row 19','row 20', ...
             'row 21','row 22','row 23','row 24','row 25','row 26','row 27','row 28','row 29','row 30', ...
             'row 31','row 32','row 33'};
columnLabels = {'maturity','horizon','rmse_RW','rmse_model'};


matrix2latex([MAT_list,HOR_list,yN_is_RW_rmse(:)*100,yN_is_rmse(:)*100], '.\results\isyN_DKW_nominal_9012.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Various Out-of-Sample Forecasts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startdaten=datenum(1990,1,1);
enddaten=datenum(2015,12,31);
data_FRBA; 
[Y,M,D]=datevec(mydate);

[logL,logL_vec,state]= logL_DKW_nominal(pvar,paras_vec,pidx,Nfac,Dt,y_vec,MATgrid,p_vec,notmissing_p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out-of-Sample Nominal Yields Forecasting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_outsample=find(Y>is_year_end);
state_oos=state(id_outsample,:);

MAT_oos=[0.5 1:10]; % bond maturity in years (6m/2y/10y)
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
[MAT_list,HOR_list,yN_oos_RW_rmse(:)*100,yN_oos_rmse(:)*100]

rowLabels = {'row 1','row 2','row 3','row 4','row 5','row 6','row 7','row 8','row 9','row 10', ...
             'row 11','row 12','row 13','row 14','row 15','row 16','row 17','row 18','row 19','row 20', ...
             'row 21','row 22','row 23','row 24','row 25','row 26','row 27','row 28','row 29','row 30', ...
             'row 31','row 32','row 33'};
columnLabels = {'maturity','horizon','rmse_RW','rmse_model'};


matrix2latex([MAT_list,HOR_list,yN_oos_RW_rmse(:)*100,yN_oos_rmse(:)*100], '.\results\oosyN_DKW_nominal_9012.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');

