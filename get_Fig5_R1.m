%% Fig 5: Ang II dose response (subcutaneous and intravenous)
%% Subcutaneous data

%%% ---- 80 ng/min (Shao et al. 2009, Shaoe et al. 2010)

AngII_circ_con_80 = 44.1;
AngII_circ_exo = 181.24;
AngII_circ_endo = 104.46;
AngII_circ_exo_SD = 35.86;
AngII_circ_endo_SD = 24.3;

fold_AngII_circ_exo = AngII_circ_exo/AngII_circ_con_80;
fold_AngII_circ_endo = AngII_circ_endo/AngII_circ_con_80;
fold_AngII_circ_80 = (AngII_circ_exo + AngII_circ_endo)/AngII_circ_con_80;
fold_AngII_circ_exo_SD = AngII_circ_exo_SD/AngII_circ_con_80;
fold_AngII_circ_endo_SD = AngII_circ_endo_SD/AngII_circ_con_80;
fold_AngII_circ_80_SD =  (AngII_circ_exo_SD + AngII_circ_endo_SD)/AngII_circ_con_80;

%%% ---- 40 ng/min (Zou et al., 1996)

fold_AngII_circ_40 = 2.936;
fold_AngII_circ_40_SD = 3.383 - 2.936;

%%% ---- 200 ng/kg/min (Campbell, 2013)

AngII_circ_con_200 = 18.3;
AngII_circ_200 = 27.3;
AngII_circ_200_SD = 5.7;
fold_AngII_circ_200 = AngII_circ_200/AngII_circ_con_200;
fold_AngII_circ_200_SD = AngII_circ_200_SD/AngII_circ_con_200;

%%% ---- 350 ng/kg/min (Campbell, 2013)

AngII_circ_con_350 = 18.3;
AngII_circ_350 = 79.5;
AngII_circ_350_SD = 28.9;
fold_AngII_circ_350 = AngII_circ_350/AngII_circ_con_350;
fold_AngII_circ_350_SD = AngII_circ_350_SD/AngII_circ_con_350;

%%% ---- 500 ng/kg/min (Campbell, 2013)

AngII_circ_con_500 = 18.3;
AngII_circ_500 = 101.9;
AngII_circ_500_SD = 15.5;
fold_AngII_circ_500 = AngII_circ_500/AngII_circ_con_500;
fold_AngII_circ_500_SD = AngII_circ_500_SD/AngII_circ_con_500;

%% Intravenous data (Pawlowski et al., 1990)

% Doses: 10, 30, 60 ng/min
AngII_circ_IV = [3.45,9.309,22.13];
AngII_circ_IV_SD = [4.229,13.03,29.221] - AngII_circ_IV;

%% Subcutaneous simulations


doses = [(0:0.1:0.9),(1:1:9),(10:10:90),(100:50:500)]; 

AngII_circ_SC = zeros(length(doses),2);

AngII_circ_ind = 3;

x0 = load('model_SS.mat').SSdata;
xp0 = zeros(length(x0),1);

for i = 1:length(doses)
    D = doses(i)
    [S,t] = run_model(13,x0,xp0,D,'ng/kg/min','SC',false);
    [~,day_7] = min(abs(t - 7*ones(size(t))));
    S_day7  = S(:,day_7);
    S_day13 = S(:,end);
    
    AngII_circ_SC(i,:) = [S_day7(AngII_circ_ind),...
                         S_day13(AngII_circ_ind)];
end

% Figure inset 
x0_sep = load('model_SS.mat').SSdata_separate;
xp0_sep = zeros(length(x0_sep), 1);
[S80,~] = run_model(13,x0_sep,xp0_sep,80,'ng/min','SC',true);

fold_AngII_circ_endo_80_sim = S80(3,end)/S80(59,1);
fold_AngII_circ_exo_80_sim = S80(42,end)/S80(59,1);
fold_AngII_circ_tot_80_sim = S80(59,end)/S80(59,1);

%% Intravenous simulations

doses_IV = [0:1:10,15:5:60];
AngII_circ_IV_sim = zeros(length(doses_IV),1);

for i = 1:length(doses_IV)
    D = doses_IV(i)
    [S,t] = run_model(1,x0,xp0,D,'ng/min','IV',false);
    [~,t_ind] = min(abs(24*t - 0.5*ones(size(t)))); 
    AngII_circ_IV_sim(i) = S(AngII_circ_ind,t_ind);
end

%% Figure

names = {'Endogenous','Exogenous','Total'};
c = summer(6);
c2 = summer(4);

figure(5)
subplot(1,2,1)
b1 = plot(doses,AngII_circ_SC(:,1)/AngII_circ_SC(1,1),'color',c2(3,:),'linewidth',1.5);
hold on
b2 = plot(doses,AngII_circ_SC(:,2)/AngII_circ_SC(1,2),'color',c2(2,:),'linewidth',1.5);
errorbar(40/0.284,fold_AngII_circ_40,fold_AngII_circ_40_SD,'d','markerfacecolor',c2(2,:),'color','k','markersize',10,'linewidth',0.8);
b4 = errorbar(80/0.284,fold_AngII_circ_80,fold_AngII_circ_80_SD,'o','color','k','markerfacecolor',c2(2,:),'markersize',10,'linewidth',0.8);
b3 = errorbar(200,fold_AngII_circ_200,fold_AngII_circ_200_SD,'o','color','k','markerfacecolor',c2(3,:),'markersize',10,'linewidth',0.8);
errorbar(350,fold_AngII_circ_350,fold_AngII_circ_350_SD,'o','color','k','markerfacecolor',c2(3,:),'markersize',10,'linewidth',0.8);
errorbar(500,fold_AngII_circ_500,fold_AngII_circ_500_SD,'o','color','k','markerfacecolor',c2(3,:),'markersize',10,'linewidth',0.8);
hold off
set(gca,'yscale','log','fontsize',14);
ylabel('[Ang II]_{circ} (log ratio to control)');
xlabel('Dose (ng/kg/min)');
legend('Day 7','Day 13');
legend([b1,b2,b3,b4],{'Simulation; day 7','Simulation; day 13', ...
    'Data; day 7','Data; day 13'},'location','southeast');
set(legend,'fontsize',14);
xlim([0,500]);
ylim([0,25]);
ttl = title('a','fontsize',20,'fontweight','bold');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';  
dim = [.3 .5 .035 .15];
annotation('rectangle',dim)
hold off

subplot(1,2,2)
plot(doses_IV,AngII_circ_IV_sim(:,1)/AngII_circ_IV_sim(1,1),'color','k','linewidth',1.5);
hold on
errorbar([10,60],[AngII_circ_IV(1),AngII_circ_IV(end)],[AngII_circ_IV_SD(1),AngII_circ_IV_SD(end)],'o','color','k',...
    'linewidth',0.8,'markerfacecolor',c(5,:),'markersize',10);
errorbar(30,AngII_circ_IV(2),AngII_circ_IV_SD(2),'d','color','k',...
    'linewidth',0.8,'markerfacecolor',c(5,:),'markersize',10);
hold off
xlabel('Dose (ng/min)');
ylabel('[Ang II]_{circ} (ratio to control)');
set(gca,'fontsize',14);
legend({'Simulation','Data'},'location','northwest');
set(legend,'fontsize',14);
ttl = title('b','fontsize',20,'fontweight','bold');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  

axes('Position',[0.15,0.7,0.15,0.2]);
hold on
box on
bar(1,fold_AngII_circ_endo_80_sim,'facealpha',0.4,'facecolor',c(6,:));
bar(1,fold_AngII_circ_endo,'barwidth',0.5,'facecolor',c(6,:));
errorbar(1,fold_AngII_circ_endo,fold_AngII_circ_endo_SD,'linestyle','none','color','k');
bar(2,fold_AngII_circ_exo_80_sim,'facealpha',0.4,'facecolor',c(5,:));
bar(2,fold_AngII_circ_exo,'barwidth',0.5,'facecolor',c(5,:));
errorbar(2,fold_AngII_circ_exo,fold_AngII_circ_exo_SD,'linestyle','none','color','k');
set(gca,'xtick',1:3,'xticklabel',names,'fontsize',12)
bar(3,fold_AngII_circ_tot_80_sim,'facealpha',0.4,'facecolor',c2(2,:));
bar(3,fold_AngII_circ_80,'barwidth',0.5,'facecolor',c2(2,:));
errorbar(3,fold_AngII_circ_80,fold_AngII_circ_80_SD,'linestyle','none','color','k');
hold off
ylim([0,max(fold_AngII_circ_80 + fold_AngII_circ_80_SD)+0.5]);
xtickangle(45)