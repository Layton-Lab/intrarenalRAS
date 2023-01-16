%% Fig 7: intrarenal Ang II validation
%% Validation data (time series)

AngII_T_con = [164,132,133,138,161];
AngII_T = [164,148,165,209,350];
AngII_T_SD = [23,15,8,10,62];
fold_AngII_T = AngII_T./AngII_T_con;
fold_AngII_T_SD = AngII_T_SD./AngII_T_con;

%% Validation simulations
% Hypothesis 1

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

% Hypothesis 2

[S402,t402] = run_model(13,x0,40,'ng/min','SC',true, {'hypothesis2'});
varargin_noRenalfb = {'all_renal'};

fold_AngII_T_sim_402 = S402(AngII_T_ind,:)/S402(AngII_T_ind,1);
fold_AngII_T_endo_sim_402 = S402(AngII_T_endo_ind,:)/S402(AngII_T_ind,1);
fold_AngII_T_exo_sim_402 = S402(AngII_T_exo_ind,:)/S402(AngII_T_ind,1);


%% Plot

c = summer(6);

figure(6)
subplot(1,2,[1,2])
hold on
plot(t40,fold_AngII_T_sim_40,'linewidth',1.5,'color','k');
plot(t40,fold_AngII_T_endo_sim_40,'linewidth',1.5,'color',c(1,:));
plot(t40,fold_AngII_T_exo_sim_40,'linewidth',1.5,'color',c(5,:));
plot(t402,fold_AngII_T_sim_402,'linewidth',1.5,'color','k','linestyle',':');
plot(t402,fold_AngII_T_endo_sim_402,'linewidth',1.5,'color',c(1,:),'linestyle',':');
plot(t402,fold_AngII_T_exo_sim_402,'linewidth',1.5,'color',c(5,:),'linestyle',':');
plot(t40_noRenalfb,fold_AngII_T_sim_noRenalfb,'--','linewidth',1.5,'color','k');
plot(t40_noRenalfb,fold_AngII_T_endo_sim_noRenalfb,'--','linewidth',1.5,'color',c(1,:));
plot(t40_noRenalfb,fold_AngII_T_exo_sim_noRenalfb,'--','linewidth',1.5,'color',c(5,:));
errorbar([0,3,7,10,13],fold_AngII_T,fold_AngII_T_SD,'d','color','k',...
    'linewidth',0.8,'markerfacecolor',c(3,:),'markersize',10);
hold off
xlabel('Time (days)');
ylabel('[Ang II]_T (ratio to control)');
set(gca,'fontsize',14);
set(legend,'fontsize',14);
legend({'[Ang II]_T','[AngII-p]_T',...
       '[AngII-i]_T','[Ang II]_T ','[AngII-p]_T ',...
     '[AngII-i]_T', '[Ang II]_T ',...
      '[AngII-p]_T','[AngII-i]_T', ...
     '[Ang II]_T (Data)'},'location','eastoutside', 'NumColumns',4, 'box','on');

box('on')

xlim([0,13]);

