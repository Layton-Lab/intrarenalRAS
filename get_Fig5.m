%% Fig 5: intrarenal Ang II fitting and validation
%% Fitting Data

% Total, exogenous, endogenous (Shao et al. 2009, Shao et al. 2010)
AngII_T_con = 206.2;
AngII_T_exo = 387.1;
AngII_T_endo = 352;
AngII_T_exo_SD = 46.2;
AngII_T_endo_SD = 76.6;

fold_AngII_T_exo = AngII_T_exo/AngII_T_con;
fold_AngII_T_endo = AngII_T_endo/AngII_T_con;
fold_AngII_T_80 = fold_AngII_T_exo + fold_AngII_T_endo;
fold_AngII_T_exo_SD = AngII_T_exo_SD/AngII_T_con;
fold_AngII_T_endo_SD = AngII_T_endo_SD/AngII_T_con;
fold_AngII_T_80_SD = fold_AngII_T_exo_SD + fold_AngII_T_endo_SD;

% Interstitial (Nishiyama et al., 2003)
AngII_Isf_con = 2.86;
AngII_Isf = 5.74;
AngII_Isf_SD = 0.26;
fold_AngII_Isf = AngII_Isf/AngII_Isf_con;
fold_AngII_Isf_SD = AngII_Isf_SD/AngII_Isf_con;

% Apical (Zhuo et al., 2002)
AngII_api_con = 37;
AngII_api = 88;
AngII_api_SD = 22;
fold_AngII_api = AngII_api/AngII_api_con;
fold_AngII_api_SD = AngII_api_SD/AngII_api_con;

% Basolateral (Zhuo et al., 2002)
AngII_baso_con = 71;
AngII_baso = 1100;
AngII_baso_SD = 283;
fold_AngII_baso = AngII_baso/AngII_baso_con;
fold_AngII_baso_SD = AngII_baso_SD/AngII_baso_con;

%% Validation data (time series)

AngII_T_con = [164,132,133,138,161];
AngII_T = [164,148,165,209,350];
AngII_T_SD = [23,15,8,10,62];
fold_AngII_T = AngII_T./AngII_T_con;
fold_AngII_T_SD = AngII_T_SD./AngII_T_con;

%% Fitting simulations

x0 = load('model_SS.mat').SSdata_separate;

[S80,t80] = run_model(13,x0,80,'ng/min','SC',true);

% Indices of key variables
AngII_T_exo_ind = 57;
AngII_T_endo_ind = 35;
AngII_Isf_ind = 64;
AT1R_AngII_cell_Tb_ind = 70;
AngII_cell_Tb_ind = 71;
AT1R_AngII_cell_Pt_ind = 66;
AngII_cell_Pt_ind = 67;

% Whole kidney
fold_AngII_T_exo_sim = S80(AngII_T_exo_ind,end)/ S80(AngII_T_endo_ind,1);
fold_AngII_T_endo_sim = S80(AngII_T_endo_ind,end)/ S80(AngII_T_endo_ind,1);
fold_AngII_T_80_sim = fold_AngII_T_exo_sim + fold_AngII_T_endo_sim;

% Interstitial
fold_AngII_Isf_sim = S80(AngII_Isf_ind,end)/S80(AngII_Isf_ind,1);

% Apical
AngII_api_con_sim = S80(AT1R_AngII_cell_Tb_ind,1) + S80(AngII_cell_Tb_ind,1);
AngII_api_sim = S80(AT1R_AngII_cell_Tb_ind,end) + S80(AngII_cell_Tb_ind,end);
fold_AngII_api_sim = AngII_api_sim/AngII_api_con_sim;

% Basolateral
AngII_baso_con_sim = S80(AT1R_AngII_cell_Pt_ind,1) + S80(AngII_cell_Pt_ind,1);
AngII_baso_sim = S80(AT1R_AngII_cell_Pt_ind,end) + S80(AngII_cell_Pt_ind,end);
fold_AngII_baso_sim = AngII_baso_sim/AngII_baso_con_sim;

%% Validation simulations

x0 = load('model_SS.mat').SSdata_separate;
[S40,t40] = run_model(13,x0,40,'ng/min','SC',true);
varargin_noRenalfb = {'all_renal'};
[S40_noRenalfb,t40_noRenalfb] = run_model(13,x0,40,'ng/min','SC',true,varargin_noRenalfb);

AngII_T_exo_ind = 57;
AngII_T_endo_ind = 35;
AngII_T_ind = 74;

fold_AngII_T_sim_40 = S40(AngII_T_ind,:)/S40(AngII_T_ind,1);
fold_AngII_T_endo_sim_40 = S40(AngII_T_endo_ind,:)/S40(AngII_T_ind,1);
fold_AngII_T_exo_sim_40 = S40(AngII_T_exo_ind,:)/S40(AngII_T_ind,1);
fold_AngII_T_sim_noRenalfb = S40_noRenalfb(AngII_T_ind,:)/S40(AngII_T_ind,1);
fold_AngII_T_endo_sim_noRenalfb = S40_noRenalfb(AngII_T_endo_ind,:)/S40(AngII_T_ind,1);
fold_AngII_T_exo_sim_noRenalfb = S40_noRenalfb(AngII_T_exo_ind,:)/S40(AngII_T_ind,1);

%% Figure

renal_names = {'Endogenous','Exogenous','Whole kidney','Interstitial','Apical','Basolateral'};

c = summer(6);

figure(5)
subplot(1,3,1)
bar(1,fold_AngII_T_endo_sim,'facecolor',c(1,:),'facealpha',0.4,'linewidth',0.8);
hold on
bar(2,fold_AngII_T_exo_sim,'facecolor',c(2,:),'facealpha',0.4,'linewidth',0.8);
bar(3,fold_AngII_T_80_sim,'facecolor',c(3,:),'facealpha',0.4,'linewidth',0.8);
bar(1,fold_AngII_T_endo,'facecolor',c(1,:),'Barwidth',0.5,'linewidth',0.8);
bar(2,fold_AngII_T_exo,'facecolor',c(2,:),'Barwidth',0.5,'linewidth',0.8);
bar(3,fold_AngII_T_80,'facecolor',c(3,:),'Barwidth',0.5,'linewidth',0.8);
errorbar([1,2,3],[fold_AngII_T_endo,fold_AngII_T_exo,fold_AngII_T_80],...
          [fold_AngII_T_endo_SD,fold_AngII_T_exo_SD,fold_AngII_T_80_SD],...
          'linestyle','none','color','k','linewidth',0.8);
bar(4,fold_AngII_Isf_sim,'facecolor',c(4,:),'facealpha',0.4,'linewidth',0.8);
bar(4,fold_AngII_Isf,'facecolor',c(4,:),'Barwidth',0.5,'linewidth',0.8);
errorbar(4,fold_AngII_Isf,fold_AngII_Isf_SD,'linestyle','none','color','k','linewidth',0.8);
bar(5,fold_AngII_api_sim,'facecolor',c(5,:),'facealpha',0.4,'linewidth',0.8);
bar(5,fold_AngII_api,'facecolor',c(5,:),'Barwidth',0.5,'linewidth',0.8);
bar(6,fold_AngII_baso_sim,'facecolor',c(6,:),'facealpha',0.4,'linewidth',0.8);
bar(6,fold_AngII_baso,'facecolor',c(6,:),'Barwidth',0.5,'linewidth',0.8);
errorbar(5,fold_AngII_api,fold_AngII_api_SD,'marker','none','linestyle','none','color','k','linewidth',0.8);
errorbar(6,fold_AngII_baso,fold_AngII_baso_SD,'marker','none','linestyle','none','color','k','linewidth',0.8);
hold off
set(gca,'xtick',1:6,'xticklabel',renal_names,'fontsize',14,'yscale','log')
ylabel('[Ang II] (log ratio to control)');
xtickangle(45)
ylim([0,21]);
ttl = title('a','fontsize',20,'fontweight','bold');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  


c = summer(6);

subplot(1,3,[2,3])
errorbar([0,3,7,10,13],fold_AngII_T,fold_AngII_T_SD,'d','color','k',...
    'linewidth',0.8,'markerfacecolor',c(3,:),'markersize',10);
hold on
plot(t40,fold_AngII_T_sim_40,'linewidth',1.5,'color','k');
plot(t40,fold_AngII_T_endo_sim_40,'linewidth',1.5,'color',c(1,:));
plot(t40,fold_AngII_T_exo_sim_40,'linewidth',1.5,'color',c(5,:));
plot(t40_noRenalfb,fold_AngII_T_sim_noRenalfb,'--','linewidth',1.5,'color','k');
plot(t40_noRenalfb,fold_AngII_T_endo_sim_noRenalfb,'--','linewidth',1.5,'color',c(1,:));
plot(t40_noRenalfb,fold_AngII_T_exo_sim_noRenalfb,'--','linewidth',1.5,'color',c(5,:));
hold off
xlabel('Time (days)');
ylabel('[Ang II]_T (ratio to control)');
set(gca,'fontsize',14);
set(legend,'fontsize',14);
legend({'[Ang II]_T (Data)','[Ang II]_T (Simulation)','[AngII-p]_T (Simulation)',...
       '[AngII-i]_T (Simulation)','[Ang II]_T (Simulation; renal fb = 0)',...
       '[AngII-p]_T (Simulation; renal fb = 0)','[AngII-i]_T (Simulation; renal fb = 0)'},'location','eastoutside');
xlim([0,13]);
ttl = title('b','fontsize',20,'fontweight','bold');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  


