close all; clear; clc;
startdaten=datenum(1990,1,1);
enddaten=datenum(2015,12,31);

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

  
% initial guess and bounds
pvarLU=[      0.8604           1.e-4  10;   % KAPPA(1,1) 
              0.1316          1.e-4  10;   % KAPPA(2,2)
              1.4528          1.e-4  10;   % KAPPA(3,3)
              -0.7358*1e-2         -100  100;   % SIGMA(2,1)
              -4.4381*1e-2         -100  100;   % SIGMA(3,1)
              -0.9482*1e-2         -100  100;   % SIGMA(3,2)
              0.0472          -1     1;   % rho0
              3.6489           -10    10;  % rho1(1)
              0.8774            -10    10;  % rho1(2)
              0.7176           -10    10;  % rho1(3)
              0.3208             -40    40;  % lambda0(1) 
              -0.4270             -40    40;  % lambda0(2)
              -1.2439             -40    40;  % lambda0(3)    
              -0.6947             -10    10;  % SIGMAlambda1(1,1)
              2.1261             -10    10;  % SIGMAlambda1(2,1)
              3.0772             -10    10;  % SIGMAlambda1(3,1)
              0.0329             -10    10;  % SIGMAlambda1(1,2)
              -0.1468             -10    10;  % SIGMAlambda1(2,2)
              -0.3630             -10    10;  % SIGMAlambda1(3,2)
              -0.0788             -10    10;  % SIGMAlambda1(1,3)
              0.5983             -10    10;  % SIGMAlambda1(2,3)
              0.1524             -10    10;  % SIGMAlambda1(3,3)         
              0.1314*1e-2        -4.e-0  4.e-0 % delta_y(1)
              0.0187*1e-2        -4.e-0  4.e-0 % delta_y(2)
              0.0654*1e-2        -4.e-0  4.e-0 % delta_y(3)
              0        -4.e-0  4.e-0 % delta_y(4)
              0.0395*1e-2        -4.e-0  4.e-0 % delta_y(5)
              0        -4.e-0  4.e-0 % delta_y(6)
              0.0531*1e-2        -4.e-0  4.e-0 % delta_y(7)
              0.1900*1e-2        -4.e-0  4.e-0 % delta_bcf(1)
              0.2979*1e-2        -4.e-0  4.e-0 % delta_bcf(2)
              0.0269             -2     2;    % rho0_pi
              -0.0209             -20     20;  % rho1_pi(1)
               0.1072             -20     20;  % rho1_pi(2)
               -0.2214           -20     20;  % rho1_pi(3)
              -0.000712             -100  100; % sigq(1)
              -0.000038             -100  100; % sigq(2)
              -0.000601             -100  100; % sigq(3)
              0.009136              -100  100; % sigqx
              0.6006*1e-2        -5   5 % delta_tips(1)
              0.4707*1e-2        -5   5 % delta_tips(2)
              0.4112*1e-2        -5   5 % delta_tips(3)
              0.2        -5   5 % delta_option(1)
              0.2        -5   5 % delta_option(2)
                 ];                 
pvar=pvarLU(:,1);
pvarL=pvarLU(:,2);
pvarU=pvarLU(:,3);
  
%load '.\results\pvar_DKW_option_9015.mat';

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


[logL,logL_vec,state,model_IE_options]= logL_DKW_option(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);

diary est_DKW_option_output.dat
  format long
  disp(pvar);
  format long
   logL
diary off

jcount=0;
while jcount < 1000
    jcount=jcount+1;

    [pvar,logL,exitflag]=fminsearch(@logL_DKW_option,pvar,OPTIONS,paras_vec,pidx,Nfac,Dt, ...              
        y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
        notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
        IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);
        
    save pvar_DKW_option pvar;

    diary est_DKW_option_output.dat
      format long
      disp(pvar);
      format long
        logL
        exitflag
      disp('-----')
    diary off

%     [logL,logL_vec,state,model_IE_options]= logL_FRBAoption(pvar,paras_vec,pidx,Nfac,Dt, ...              
%         y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
%         notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
%         IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);
%     
%     id_options=find(notmissing_options); 
%     IE_options_avg=[mean(IE_options(:,1:3)')',mean(IE_options(:,4:6)')'];
%     subplot(2,1,1); plot([model_IE_options(id_options,1),IE_options_avg(id_options,1)]);
%     subplot(2,1,2); plot([model_IE_options(id_options,2),IE_options_avg(id_options,2)]);

end  

diary est_DKW_option_output.dat
 disp('****************')
diary off

[logL,logL_vec,state,model_IE_options]= logL_DKW_option(pvar,paras_vec,pidx,Nfac,Dt, ...              
    y_vec,MATgrid,ytips_vec,TIPSgrid,p_vec,bcf_vec,bcfLT_vec,horLT, ...
    notmissing_tips,notmissing_p,notmissing_bcf,notmissing_bcfLT, ...
    IE_options,MATgrid_options,STRIKEgrid_options,notmissing_options);
