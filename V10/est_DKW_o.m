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

% load '.\results\pvar_DKW_9012.mat';
% pvar=[pvar;0.005;0.005];

%load '.\results\pvar_DKW_o_9012.mat';
pvar = [
   0.879827174869795
   0.136676781421871
   1.447294240878790
  -0.006935343022228
  -0.045186863363826
  -0.009686920348521
   0.047345903856018
   3.772343300716693
   0.895183130460821
   0.739790198129056
   0.307279490317847
  -0.408537703642677
  -1.233195040745562
  -0.694010465657616
   2.113043355724217
   3.178832667409866
   0.030392958476213
  -0.143960242446123
  -0.350526171296919
  -0.066405899989179
   0.570349172643620
   0.128095773976990
   0.001323657015903
   0.000183208888583
   0.000657008625219
   0.000000000000000
   0.000397218112930
  -0.000000000000001
   0.000529298428788
   0.001886479136480
   0.002944078095088
   0.026807560337462
  -0.098763320389107
   0.118172707006635
  -0.229943686982747
  -0.000657406648920
  -0.000055579646865
  -0.000615032637428
   0.009139346059945
   0.005986509410064
   0.004702858981277
   0.004114711779935
   0.005000000000000
   0.005000000000000]; 

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


itersize=1000;
Tol=1.e-9;
Tol2=1.e-9;

OPTIONS = optimset('Display','iter','Diagnostics','off',...
                   'MaxFunEvals',itersize);
OPTIONS = optimset(OPTIONS,'TolFun', Tol, 'TolX',Tol2);
OPTIONS = optimset(OPTIONS,'MaxSQPIter',100);


[logL,logL_vec,state,model_IE_options]= logL_DKW_o(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);

diary est_DKW_o_9012.dat
  format long
  disp(pvar);
  format long
   logL
diary off

jcount=0;
while jcount < 1000
    jcount=jcount+1;

    [pvar,logL,exitflag]=fminsearch(@logL_DKW_o,pvar,OPTIONS,paras_vec,pidx,Nfac,Dt, ...              
        y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
        notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
        IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);
        
    save pvar_DKW_o_9012 pvar;

    diary est_DKW_o_9012.dat
      format long
      disp(pvar);
      format long
        logL
        exitflag
      disp('-----')
    diary off

end  

diary est_DKW_o_9012.dat
 disp('****************')
diary off

[logL,logL_vec,state,model_IE_options]= logL_DKW_o(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);
