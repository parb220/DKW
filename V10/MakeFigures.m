close all; clear; clc;
startdaten=datenum(1990,1,1);
enddaten=datenum(2015,12,31);
data_FRBA; 
[Y,M,D]=datevec(mydate); 
dateaxis=Y+((M-1)*30+D)/365;

% % % Figure: NOMINAL Yield Curves %
% % fs=14;
% % load ('..\data\allyield.txt');
% % allyield(allyield(:,1)<1990,:)=[];
% % yields=allyield(:,4:end)/100; % 
% % MATgrid_allyield=[0.25 0.5 1:30];      % maturities for nominal yields
% % date_axis=allyield(:,1)+((allyield(:,2)-1)*30+allyield(:,3))/365;
% % fig_nomyc=figure;
% % mesh(date_axis,MATgrid_allyield,yields');
% % view([45,45])
% % title('Nominal Yield Curves: 1990-2015','FontSize',fs); 
% % xlabel('date','FontSize',fs);
% % ylabel('maturity','FontSize',fs);
% % zlabel('yield','FontSize',fs);
% % 
% % print(fig_nomyc,'-depsc2','..\..\Drafts\nomyc.eps');
% % 
% % 
% % % Figure: REAL Yield Curves %
% % load ('..\data\alltips.txt');
% % alltips(alltips(:,1)<1990,:)=[];
% % yields_alltips=alltips(:,4:end)/100;
% % MATgrid_alltips=[2:20];
% % id_beg=min(find(yields_alltips(:,4)));
% % yields_alltips(1:(id_beg-1),:)=NaN(id_beg-1,19);
% % yields_alltips(find(yields_alltips(:,1)==0),1)=NaN;
% % yields_alltips(find(yields_alltips(:,2)==0),2)=NaN;
% % yields_alltips(find(yields_alltips(:,3)==0),3)=NaN;
% % date_axis_alltips=alltips(:,1)+((alltips(:,2)-1)*30+alltips(:,3))/365;
% % fig_realyc=figure;
% % mesh(date_axis_alltips,MATgrid_alltips,yields_alltips');
% % view([45,45])
% % title('Real Yield Curves: 1999-2015','FontSize',fs); 
% % xlabel('date','FontSize',fs);
% % ylabel('maturity','FontSize',fs);
% % zlabel('yield','FontSize',fs);
% % print(fig_realyc,'-depsc2','..\..\Drafts\realyc.eps');

fs=12;
lw=3;

% % Figure: 5y-on-5y Breakeven Inflation Rate % 
fig_BEI=figure;
subplot(2,1,1);
plot(dateaxis(id_tips),[yields(id_tips,7),alltips(id_tips,4)],'LineWidth',lw);
legend('Nominal','Real');
title('5-year Treasuries and TIPS yields','FontSize',fs);

subplot(2,1,2);
plot(dateaxis(id_tips),yields(id_tips,7)-alltips(id_tips,4),'b','LineWidth',lw);
title('5-year Breakeven Inflation Rate','FontSize',fs);
print(fig_BEI,'-depsc2','..\..\Drafts\BEI.eps');

% % Figure: Option-Implied Inflation Expectations 
fig_IEoption=figure;
IE_options_avg=[mean(IE_options(:,1:3)')',mean(IE_options(:,4:6)')'];
BEI_3y=yields(:,5)-alltips(:,2);

subplot(2,1,1);
plot(dateaxis(id_options),IE_options_avg(id_options,2),'b',dateaxis(id_options),IE_options_avg(id_options,1),'b--','LineWidth',lw);
legend('3-year option-implied IE','1-year option-implied IE');
title('Option-Implied Inflation Expecations','FontSize',fs);

subplot(2,1,2);
plot(dateaxis(id_options),IE_options_avg(id_options,2),'b',dateaxis(id_options),BEI_3y(id_options),'r--','LineWidth',lw);
legend('3-year option-implied IE','3-year BEI');
title('Option-Implied Inflation Expecation vs BEI','FontSize',fs);
print(fig_IEoption,'-depsc2','..\..\Drafts\IEoption.eps');


% % Figure: state variables 
[logL_DKW,state_DKW,IE_DKW,IRP_DKW,caps_DKW,floors_DKW]=analysis_DKW;
[logL_DKW_o,state_DKW_o,IE_DKW_o,IRP_DKW_o,caps_DKW_o,floors_DKW_o]=analysis_DKW_o;
[logL_DKWv_o,state_DKWv_o,IE_DKWv_o,IRP_DKWv_o,caps_DKWv_o,floors_DKWv_o]=analysis_DKWv_o;

fig_state=figure;
subplot(2,2,1);
plot(dateaxis,state_DKWv_o(:,2),'b',dateaxis,state_DKW_o(:,2),'g--',dateaxis,state_DKW(:,2),'r.','LineWidth',lw);
legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
title('$x_1$','FontSize',fs,'Interpreter','Latex');

subplot(2,2,2);
plot(dateaxis,state_DKWv_o(:,3),'b',dateaxis,state_DKW_o(:,3),'g--',dateaxis,state_DKW(:,3),'r.','LineWidth',lw);
legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
title('$x_2$','FontSize',fs,'Interpreter','Latex');

subplot(2,2,3);
plot(dateaxis,state_DKWv_o(:,4),'b',dateaxis,state_DKW_o(:,4),'g--',dateaxis,state_DKW(:,4),'r.','LineWidth',lw);
legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
title('$x_3$','FontSize',fs,'Interpreter','Latex');

subplot(2,2,4);
plot(dateaxis,state_DKWv_o(:,5),'b','LineWidth',lw);
legend('Model^{NRO}_{xv}');
title('$v$','FontSize',fs,'Interpreter','Latex');
print(fig_state,'-depsc2','..\..\Drafts\state.eps');

% FIGURE: cap prices %
fig_caps=figure;
subplot(2,3,1); 
plot(dateaxis(id_caps),exp(caps_vec(id_caps,1)),'k',dateaxis(id_caps),caps_DKWv_o(id_caps,1),'b',dateaxis(id_caps),caps_DKW_o(id_caps,1),'g--',dateaxis(id_caps),caps_DKW(id_caps,1),'r.','LineWidth',lw);
title('1-year at 1% strike','FontSize',fs);
% xlim([dateaxis(id_caps(1)) dateaxis(id_caps(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('1-year CAPS with 1% strike','FontSize',fs);


subplot(2,3,2); 
plot(dateaxis(id_caps),exp(caps_vec(id_caps,2)),'k',dateaxis(id_caps),caps_DKWv_o(id_caps,2),'b',dateaxis(id_caps),caps_DKW_o(id_caps,2),'g--',dateaxis(id_caps),caps_DKW(id_caps,2),'r.','LineWidth',lw);
title('1-year at 2% strike','FontSize',fs);
% xlim([dateaxis(id_caps(1)) dateaxis(id_caps(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('1-year CAPS with 2% strike','FontSize',fs);

subplot(2,3,3); 
plot(dateaxis(id_caps),exp(caps_vec(id_caps,3)),'k',dateaxis(id_caps),caps_DKWv_o(id_caps,3),'b',dateaxis(id_caps),caps_DKW_o(id_caps,3),'g--',dateaxis(id_caps),caps_DKW(id_caps,3),'r.','LineWidth',lw);
title('1-year at 3% strike','FontSize',fs);
% xlim([dateaxis(id_caps(1)) dateaxis(id_caps(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('1-year CAPS with 3% strike','FontSize',fs);

subplot(2,3,4); 
plot(dateaxis(id_caps),exp(caps_vec(id_caps,4)),'k',dateaxis(id_caps),caps_DKWv_o(id_caps,4),'b',dateaxis(id_caps),caps_DKW_o(id_caps,4),'g--',dateaxis(id_caps),caps_DKW(id_caps,4),'r.','LineWidth',lw);
title('3-year at 1% strike','FontSize',fs);
% xlim([dateaxis(id_caps(1)) dateaxis(id_caps(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('3-year CAPS with 1% strike','FontSize',fs);

subplot(2,3,5); 
plot(dateaxis(id_caps),exp(caps_vec(id_caps,5)),'k',dateaxis(id_caps),caps_DKWv_o(id_caps,5),'b',dateaxis(id_caps),caps_DKW_o(id_caps,5),'g--',dateaxis(id_caps),caps_DKW(id_caps,5),'r.','LineWidth',lw);
title('3-year at 2% strike','FontSize',fs);
% xlim([dateaxis(id_caps(1)) dateaxis(id_caps(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('3-year CAPS with 2% strike','FontSize',fs);

subplot(2,3,6); 
plot(dateaxis(id_caps),exp(caps_vec(id_caps,6)),'k',dateaxis(id_caps),caps_DKWv_o(id_caps,6),'b',dateaxis(id_caps),caps_DKW_o(id_caps,6),'g--',dateaxis(id_caps),caps_DKW(id_caps,6),'r.','LineWidth',lw);
title('3-year at 3% strike','FontSize',fs);
% xlim([dateaxis(id_caps(1)) dateaxis(id_caps(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('3-year CAPS with 3% strike','FontSize',fs);

print(fig_caps,'-depsc2','..\..\Drafts\caps.eps');


% FIGURE: floor prices %
fig_floors=figure;
subplot(2,3,1); 
plot(dateaxis(id_floors),exp(floors_vec(id_floors,1)),'k',dateaxis(id_floors),floors_DKWv_o(id_floors,1),'b',dateaxis(id_floors),floors_DKW_o(id_floors,1),'g--',dateaxis(id_floors),floors_DKW(id_floors,1),'r.','LineWidth',lw);
title('1-year at 1% strike','FontSize',fs);
% xlim([dateaxis(id_floors(1)) dateaxis(id_floors(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('1-year FLOORS  with 1% strike','FontSize',fs);

subplot(2,3,2); 
plot(dateaxis(id_floors),exp(floors_vec(id_floors,2)),'k',dateaxis(id_floors),floors_DKWv_o(id_floors,2),'b',dateaxis(id_floors),floors_DKW_o(id_floors,2),'g--',dateaxis(id_floors),floors_DKW(id_floors,2),'r.','LineWidth',lw);
title('1-year at 2% strike','FontSize',fs);
% xlim([dateaxis(id_floors(1)) dateaxis(id_floors(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('1-year FLOORS with 2% strike','FontSize',fs);

subplot(2,3,3); 
plot(dateaxis(id_floors),exp(floors_vec(id_floors,3)),'k',dateaxis(id_floors),floors_DKWv_o(id_floors,3),'b',dateaxis(id_floors),floors_DKW_o(id_floors,3),'g--',dateaxis(id_floors),floors_DKW(id_floors,3),'r.','LineWidth',lw);
title('1-year at 3% strike','FontSize',fs);
% xlim([dateaxis(id_floors(1)) dateaxis(id_floors(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('1-year FLOORS with 3% strike','FontSize',fs);

subplot(2,3,4); 
plot(dateaxis(id_floors),exp(floors_vec(id_floors,4)),'k',dateaxis(id_floors),floors_DKWv_o(id_floors,4),'b',dateaxis(id_floors),floors_DKW_o(id_floors,4),'g--',dateaxis(id_floors),floors_DKW(id_floors,4),'r.','LineWidth',lw);
title('3-year at 1% strike','FontSize',fs);
% xlim([dateaxis(id_floors(1)) dateaxis(id_floors(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('3-year FLOORS with 1% strike','FontSize',fs);

subplot(2,3,5); 
plot(dateaxis(id_floors),exp(floors_vec(id_floors,5)),'k',dateaxis(id_floors),floors_DKWv_o(id_floors,5),'b',dateaxis(id_floors),floors_DKW_o(id_floors,5),'g--',dateaxis(id_floors),floors_DKW(id_floors,5),'r.','LineWidth',lw);
title('3-year at 2% strike','FontSize',fs);
% xlim([dateaxis(id_floors(1)) dateaxis(id_floors(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('3-year FLOORS with 2% strike','FontSize',fs);

subplot(2,3,6); 
plot(dateaxis(id_floors),exp(floors_vec(id_floors,6)),'k',dateaxis(id_floors),floors_DKWv_o(id_floors,6),'b',dateaxis(id_floors),floors_DKW_o(id_floors,6),'g--',dateaxis(id_floors),floors_DKW(id_floors,6),'r.','LineWidth',lw);
title('3-year at 3% strike','FontSize',fs);
% xlim([dateaxis(id_floors(1)) dateaxis(id_floors(end))]);
% legend('Data','Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
% title('3-year FLOORS with 3% strike','FontSize',fs);

print(fig_floors,'-depsc2','..\..\Drafts\floors.eps');


fig_alldecomp=figure;
subplot(2,2,1); plot(dateaxis,IE_DKWv_o(:,1),'b',dateaxis,IE_DKW_o(:,1),'g--',dateaxis,IE_DKW(:,1),'r.','LineWidth',lw); 
title('1-year Inflation Expectations','FontSize',fs); 
% legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
xlim([dateaxis(1) dateaxis(end)])

subplot(2,2,3); plot(dateaxis,IE_DKWv_o(:,2),'b',dateaxis,IE_DKW_o(:,2),'g--',dateaxis,IE_DKW(:,2),'r.','LineWidth',lw); 
title('5-year Inflation Expectations','FontSize',fs); 
% legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
xlim([dateaxis(1) dateaxis(end)])

subplot(2,2,2); plot(dateaxis,IRP_DKWv_o(:,1),'b',dateaxis,IRP_DKW_o(:,1),'g--',dateaxis,IRP_DKW(:,1),'r.','LineWidth',lw); 
title('1-year Inflation Risk Premiums','FontSize',fs); 
% legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
xlim([dateaxis(1) dateaxis(end)])

subplot(2,2,4); plot(dateaxis,IRP_DKWv_o(:,2),'b',dateaxis,IRP_DKW_o(:,2),'g--',dateaxis,IRP_DKW(:,2),'r.','LineWidth',lw); 
title('5-year Inflation Risk Premiums','FontSize',fs); 
% legend('Model^{NRO}_{xv}','Model^{NRO}_{x}','Model^{NR}_{x}');
xlim([dateaxis(1) dateaxis(end)])

print(fig_alldecomp,'-depsc2','..\..\Drafts\alldecomp.eps');


fig_allIE=figure;
load '..\data\umich.txt';
id_mich=find(umich(:,4)>0 & datenum(umich(:,1:3))>mydate(1) & datenum(umich(:,1:3))<=mydate(end));

load '..\data\FMAieS.txt';
id_ieS=find(FMAieS(:,4)>0 & datenum(FMAieS(:,1:3))>mydate(1) & datenum(FMAieS(:,1:3))<=mydate(end));

plot(dateaxis,IE_DKWv_o(:,1),'b-','LineWidth',lw); hold on;
plot(umich(id_mich,1)+umich(id_mich,2)/12+umich(id_mich,3)/365,umich(id_mich,4)/100,'r-'); hold on;
plot(FMAieS(id_ieS,1)+FMAieS(id_ieS,2)/12+FMAieS(id_ieS,3)/365,FMAieS(id_ieS,4)/100,'ko'); hold off;
legend('Model^{NRO}_{xv}','Michigan Survey','Blue Chip Survey');
title('1-year Inflation Expectation','FontSize',fs);
xlim([dateaxis(1) dateaxis(end)]);
print(fig_allIE,'-depsc2','..\..\Drafts\allIE_FRBAfull.eps');

