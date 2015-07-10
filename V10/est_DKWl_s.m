close all; clear; clc;
startdaten=datenum(1990,1,1);
enddaten=datenum(2012,12,31);

data_FRBA; 

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

SIGMA_L=0.01;

theta=zeros(Nfac,1);

delta_bcfLT=0.0075;

delta_p=0;

    
paras_vec=[reshape(KAPPA,Nfac^2,1);reshape(SIGMA,Nfac^2,1);theta;
      rho0;rho1;lambda0;reshape(SIGMAlambda1,Nfac^2,1);
      delta_y;delta_bcf;delta_bcfLT;rho0_pi;rho1_pi;sigq;
      sigqx;delta_p;delta_tips;delta_swap;KAPPA_L;SIGMA_L;theta_L;
      rho1_L;rhoL_L;lambda0_L;SIGMAlambda1_L];

load '.\results\pvar_DKWl_9012.mat';
pvar=[pvar(1:42);0.01;0.01;0.01;pvar(43:end)];

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


itersize=3000;
Tol=1.e-9;
Tol2=1.e-9;

OPTIONS = optimset('Display','iter','Diagnostics','off',...
                   'MaxFunEvals',itersize);
OPTIONS = optimset(OPTIONS,'TolFun', Tol, 'TolX',Tol2);
OPTIONS = optimset(OPTIONS,'MaxSQPIter',100);


[logL,logL_vec,state]= logL_DKWl_s(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,yS_vec,notmissing_swap);

diary est_DKWl_s.dat
  format long
   disp(pvar)
  format long
   logL
diary off

jcount=0;
while jcount < 1000
  jcount=jcount+1;

  [pvar,logL,exitflag]=fminsearch('logL_DKWl_s',pvar,OPTIONS,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,yS_vec,notmissing_swap);

  save pvar_DKWl_s pvar;
  
  diary est_DKWl_s.dat
      format long
        disp(pvar)
      format long
        logL
        exitflag
      disp('-----')
  diary off


end  

diary est_DKWl_s.dat
 disp('****************')
diary off

[logL,logL_vec,state]= logL_DKWl_s(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT,yS_vec,notmissing_swap);




