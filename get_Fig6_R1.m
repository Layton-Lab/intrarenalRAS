%% Fig 6: intrarenal Ang II fitting
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

%% Fitting simulations

x0 = load('model_SS.mat').SSdata_separate;
[S80,t80] = run_model(13,x0,80,'ng/min','SC',true);

x0 = load('model_SS.mat').SSdata_separate;
[S802,t802] = run_model(13,x0,80,'ng/min','SC',true,{'hypothesis2'});

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
fold_AngII_T_exo_sim2 = S802(AngII_T_exo_ind,end)/ S802(AngII_T_endo_ind,1);
fold_AngII_T_endo_sim2 = S802(AngII_T_endo_ind,end)/ S802(AngII_T_endo_ind,1);
fold_AngII_T_80_sim2 = fold_AngII_T_exo_sim2 + fold_AngII_T_endo_sim2;

% Interstitial
fold_AngII_Isf_sim = S80(AngII_Isf_ind,end)/S80(AngII_Isf_ind,1);
fold_AngII_Isf_sim2 = S802(AngII_Isf_ind,end)/S802(AngII_Isf_ind,1);

% Apical
AngII_api_con_sim = S80(AT1R_AngII_cell_Tb_ind,1) + S80(AngII_cell_Tb_ind,1);
AngII_api_sim = S80(AT1R_AngII_cell_Tb_ind,end) + S80(AngII_cell_Tb_ind,end);
fold_AngII_api_sim = AngII_api_sim/AngII_api_con_sim;
AngII_api_con_sim2 = S802(AT1R_AngII_cell_Tb_ind,1) + S802(AngII_cell_Tb_ind,1);
AngII_api_sim2 = S802(AT1R_AngII_cell_Tb_ind,end) + S802(AngII_cell_Tb_ind,end);
fold_AngII_api_sim2 = AngII_api_sim2/AngII_api_con_sim2;

% Basolateral
AngII_baso_con_sim = S80(AT1R_AngII_cell_Pt_ind,1) + S80(AngII_cell_Pt_ind,1);
AngII_baso_sim = S80(AT1R_AngII_cell_Pt_ind,end) + S80(AngII_cell_Pt_ind,end);
fold_AngII_baso_sim = AngII_baso_sim/AngII_baso_con_sim;
AngII_baso_con_sim2 = S802(AT1R_AngII_cell_Pt_ind,1) + S802(AngII_cell_Pt_ind,1);
AngII_baso_sim2 = S802(AT1R_AngII_cell_Pt_ind,end) + S802(AngII_cell_Pt_ind,end);
fold_AngII_baso_sim2 = AngII_baso_sim2/AngII_baso_con_sim2;

exo = [fold_AngII_T_exo fold_AngII_T_exo_sim fold_AngII_T_exo_sim2];
endo = [fold_AngII_T_endo fold_AngII_T_endo_sim fold_AngII_T_endo_sim2];
T = [fold_AngII_T_80 fold_AngII_T_80_sim fold_AngII_T_80_sim2];
isf = [fold_AngII_Isf fold_AngII_Isf_sim fold_AngII_Isf_sim2];
api = [fold_AngII_api fold_AngII_api_sim fold_AngII_api_sim2];
baso = [fold_AngII_baso fold_AngII_baso_sim fold_AngII_baso_sim2];
y = [endo; exo; T; isf;api;baso];

%% Plot (new)

renal_names = {'Endogenous','Exogenous','Whole kidney','Interstitial','Apical','Basolateral'};

c = summer(6);
idxs = [1 1 5];
x = 0.225; 
figure(6)
b = bar(y,'FaceColor','flat');
hold on
errorbar([1-x,2-x,3-x],[fold_AngII_T_endo,fold_AngII_T_exo,fold_AngII_T_80],...
          [fold_AngII_T_endo_SD,fold_AngII_T_exo_SD,fold_AngII_T_80_SD],...
          'linestyle','none','color','k','linewidth',0.8);
errorbar(4-x,fold_AngII_Isf,fold_AngII_Isf_SD,'linestyle','none','color','k','linewidth',0.8);
errorbar(5-x,fold_AngII_api,fold_AngII_api_SD,'marker','none','linestyle','none','color','k','linewidth',0.8);
errorbar(6-x,fold_AngII_baso,fold_AngII_baso_SD,'marker','none','linestyle','none','color','k','linewidth',0.8);
hold off
set(gca,'xtick',1:6,'xticklabel',renal_names,'fontsize',14,'yscale','log')
ylabel('[Ang II] (log ratio to control)');
xtickangle(45)
ylim([0,21]);
b(1).CData = [0 0 0];
for k = 2:size(y,2)
    b(k).CData = c(idxs(k),:);
end
legend('Data','Mechanism (i)','Mechanism (ii)')
set(legend,'fontsize',12,'loc','northwest')

