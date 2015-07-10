Nfac=3;     % number of factors
Dt=7/365;  % weekly data

load('../data/FMAyield.txt'); yield=FMAyield;
load('../data/FMAcpi.txt'); cpi=FMAcpi;
load('../data/FMAtips.txt'); tips=FMAtips;
% load('../data/FMAtipsCA.txt'); tips=FMAtipsCA;
load('../data/FMAbcf.txt'); bcf=FMAbcf;
load('../data/FMAbcfLT.txt'); bcfLT=FMAbcfLT;
load('../data/caps.txt');
load('../data/floors.txt');
load('../data/oil.txt');

load ../data/FMAieS.txt; ieS=FMAieS;
load ../data/FMAieL.txt; ieL=FMAieL;

alldaten=datenum(yield(:,1:3));
id_data=find(alldaten>=startdaten & alldaten<=enddaten);
id_beg=min(id_data); id_end=max(id_data);

% id_data=find(yield(:,1)>=startyear & yield(:,1)<=endyear);
mydate=datenum(yield(id_beg:id_end,1:3));


MATgrid=[0.25 0.5 1 2 4 7 10];      % maturities for nominal yields
Ny=length(MATgrid);
y_vec=yield(id_beg:id_end,4:10);    % maturities in MATgrid
y_vec=y_vec/100;
T=size(y_vec,1);

TIPSgrid=[5 7 10];
Nytips=length(TIPSgrid);
ytips_vec=tips(id_beg:id_end,4:6);
ytips_vec=ytips_vec/100;
notmissing_tips=1-(ytips_vec(:,1)==0);          % 1, if not missing; 0, otherwise
id_tips=find(notmissing_tips==1);

p_vec0=cpi(id_beg:id_end,4);
notmissing_p=1-(p_vec0==0);          % 1, if not missing; 0, otherwise
id_p=find(notmissing_p==1);
p_vec=zeros(size(notmissing_p));
p_vec(id_p)=log(p_vec0(id_p));

MATgrid_caps=[1 3];
STRIKEgrid_caps=[0.01 0.02 0.03]; 
Ncaps=length(MATgrid_caps)*length(STRIKEgrid_caps);

caps_vec=caps(id_beg:id_end,4:9);   % [1yr_1%,1yr_2%,1yr_3%,3yr_1%,3yr_2%,3yr_3%]
caps_vec=caps_vec/10^4;
notmissing_caps=1-(caps_vec(:,1)==0);
id_caps=find(notmissing_caps==1);
caps_vec(id_caps,:)=log(caps_vec(id_caps,:));

MATgrid_floors=[1 3];
STRIKEgrid_floors=[0.01 0.02 0.03]; 
Nfloors=length(MATgrid_floors)*length(STRIKEgrid_floors);

floors_vec=floors(id_beg:id_end,4:9);   % [1yr_1%,1yr_2%,1yr_3%,3yr_1%,3yr_2%,3yr_3%]
floors_vec=floors_vec/10^4;
notmissing_floors=1-(floors_vec(:,1)==0);
id_floors=find(notmissing_floors==1);
floors_vec(id_floors,:)=log(floors_vec(id_floors,:));

% construct options-implied inflation expectations
id_options=max(min(id_caps),min(id_floors)):T; id_options=id_options';
MATgrid_options=[1 3];
STRIKEgrid_options=[0.01 0.02 0.03]; 
Noptions=length(MATgrid_options)*length(STRIKEgrid_options);
load ('../data/allyield.txt');
yvec_options=allyield(id_beg:id_end,[6 8])/100; % 
IE_options=zeros(T,Noptions);
for i=1:length(MATgrid_options)
    for j=1:length(STRIKEgrid_options)
        tmpID=(i-1)*length(STRIKEgrid_options)+j;
        MAT=MATgrid_options(i); STRIKE=STRIKEgrid_options(j);
        tmp_IE_options=(exp(caps_vec(id_options,tmpID))-exp(floors_vec(id_options,tmpID))) ...
                                     ./(exp(-MAT*yvec_options(id_options,i)))+(1+STRIKE)^MAT;
        IE_options(id_options,tmpID)=1/MAT*log(tmp_IE_options);                                 
    end
end
notmissing_options=1-(IE_options(:,1)==0);



MATgrid_oil=oil(id_beg:id_end,9:11)'/360;
oil_vec=oil(id_beg:id_end,5:7);   % [F2, F4, F13] contracts
notmissing_oil=1-(oil_vec(:,1)==0);
id_oil=find(notmissing_oil==1);
oil_vec(id_oil,:)=log(oil_vec(id_oil,:));

bcf05=bcf(id_beg:id_end,4);  % forecast of 5-year yield
bcf10=bcf(id_beg:id_end,5);  % forecast of 10-year yield
bcf6_11=bcfLT(id_beg:id_end,4);
bcf6_11_year=bcfLT(id_beg:id_end,5);

id_bcfLT=find(bcf6_11);
bcf6_11_date=[bcf6_11_year(id_bcfLT),12*ones(length(id_bcfLT),1),31*ones(length(id_bcfLT),1)];
hortmp= (datenum(bcf6_11_date)-mydate(id_bcfLT))/365.25;
hortmp= hortmp-5.0;
horLT= zeros(size(bcf6_11));
horLT(id_bcfLT)=hortmp;

bcf_vec= [bcf05 bcf10]/100;
bcfLT_vec=bcf6_11/100;
notmissing_bcf= 1- (bcf_vec(:,1) == 0);  % =1 if nonmissing, =0 if missing
notmissing_bcfLT= 1 - (bcfLT_vec == 0);  % =1 if nonmissing, =0 if missing


ieS_vec=ieS(id_beg:id_end,end);
ieS_vec=ieS_vec/100;
notmissing_ieS=1-(ieS_vec==0);  
id_ieS=find(notmissing_ieS==1);

ieL_vec=ieL(id_beg:id_end,end);
ieL_vec=ieL_vec/100;
notmissing_ieL=1-(ieL_vec==0);  
id_ieL=find(notmissing_ieL==1);


% % % convenience yields in oil futures %
% % load ('..\data\allyield.txt');
% % allyield(allyield(:,1)<1990,:)=[];
% % yields=allyield(:,4:end)/100; % 
% % MATgrid_allyield=[0.25 0.5 1:30];      % maturities for nominal yields
% % 
% % load ('..\data\oil.txt');
% % oil(oil(:,1)<1990,:)=[];
% % oil_vec=log(oil(:,4:7));   % [F1, F2, F4, F13] contracts
% % MATgrid_oil=oil(:,8:11)'/360;
% % 
% % T=size(yields,1);
% % nom_rate=zeros(T,1);
% % for i=1:T   
% %     nom_rate(i)=interp1(MATgrid_allyield,yields(i,:),MATgrid_oil(4,i)); 
% % end
% % 
% % cy_vec=(exp(nom_rate.*MATgrid_oil(4,:)').*exp(oil_vec(:,1))-exp(oil_vec(:,4)))./exp(oil_vec(:,1));
% % s_vec=oil_vec(:,1);
% % % rho_ds=corrcoef(cy_vec(2:end)-cy_vec(1:(end-1)),s_vec(2:end)-s_vec(1:(end-1)));
% % 
% % 
% % load('..\results\state_FRBA_9015.mat');
% % % B=[rho0;rho1;rhov;rhod;rhos]
% % [B,BINT,R]=regress(yields(:,1),[ones(T,1),state_FRBA_9015(:,5:8),cy_vec,s_vec]);
% % 
% % % B1=[phi0,phid]; phi1=phiv=phis=0 as initial values
% % [B1,BINT1,R1]=regress(cy_vec(2:end)-cy_vec(1:(end-1)),[ones(T-1,1),cy_vec(1:(end-1))]*Dt);
% % 
% % [B2,BINT2,R2]=regress(s_vec(2:end)-s_vec(1:(end-1)),[ones(T-1,1),state_FRBA_9015(1:(end-1),5:8),cy_vec(1:(end-1)),s_vec(1:(end-1))]*Dt);
% % % SIGMAs=sqrt(var(R1));

% forehor=6;
% [B,BINT,R,RINT,STATS] =regress(p_vec(id_p((1+forehor):end))-p_vec(id_p(1:(end-forehor))),[ones(length(id_p)-forehor,1),cy_vec(id_p(1:(end-forehor)))]);
% convenience yields in oil futures %

load ('../data/allyield.txt');
yields=allyield(id_beg:id_end,4:end)/100; % 
MATgrid_allyield=[0.25 0.5 1:30];      % maturities for nominal yields

load ('../data/alltips.txt');
alltips=alltips(id_beg:id_end,4:end)/100; % 
MATgrid_alltips=[2:30];    

load ('../data/swap.txt');
swap=swap(id_beg:id_end,4:end); % 
MATgrid_swap=[1:10];     
swap_vec=swap(:,TIPSgrid); % same maturities TIPSgrid=[5 7 10];
notmissing_swap= 1-(swap_vec==999);  % =1 if nonmissing, =0 if missing
Nswap=size(swap_vec,2);
yS_vec=zeros(T,Nswap);
for i=1:Nswap
    MAT=TIPSgrid(i);
    id_allyield=find(MATgrid_allyield==MAT);
    tmp_id_swap=find(notmissing_swap(:,i)==1);
    swap_vec(tmp_id_swap,i)=swap_vec(tmp_id_swap,i)/100;
    yS_vec(tmp_id_swap,i)=yields(tmp_id_swap,id_allyield)-swap_vec(tmp_id_swap,i);
end



load ('../data/oil.txt');
oil_vec2=log(oil(id_beg:id_end,4:7));   % [F1, F2, F4, F13] contracts
MATgrid_oil2=oil(id_beg:id_end,8:11)'/360;

% nom_rate=zeros(T,1);
% for i=1:T   
%     nom_rate(i)=interp1(MATgrid_allyield,yields(i,:),MATgrid_oil2(4,i)); 
% end
% 
% cy_vec=(exp(nom_rate.*MATgrid_oil2(4,:)').*exp(oil_vec2(:,1))-exp(oil_vec2(:,4)))./exp(oil_vec2(:,1));
% s_vec=oil_vec2(:,1);
