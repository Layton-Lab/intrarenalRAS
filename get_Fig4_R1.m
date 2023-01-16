%% Fig 4 - AGT, PRA, and Ang I time series
%% Data

% AngI_circ
fold_AngI_circ = [1,0.26,0.32,0.38,0.165];
fold_AngI_circ_SD = [1,0.305,0.365,0.46,0.215] - fold_AngI_circ;

% AngI_T
AngI_T_con = [130,130,128,112,122];
AngI_T = [130,121,105,107,104];
AngI_T_SD = [130,131,117,120,120]-AngI_T;
fold_AngI_T = AngI_T./AngI_T_con;
fold_AngI_T_SD = AngI_T_SD./AngI_T_con;

% PRA
fold_PRA = [1.07/4.5,0.52/4.75,0.64/4.34,0.05/4.17];
fold_PRA_SD = [0.20/4.5,0.21/4.75,0.14/4.34,0.02/4.17];

% AGT
fold_AGT = 1.58;
fold_AGT_SD = 1.75 - 1.58;

%% Simulation

x0 = load('model_SS.mat').SSdata;
varargin_noAGTfb = {'circ_AGT'};
varargin_noACEfb = {'circ_ACE'};
varargin_noReninfb = {'renin'};

[S40,t40] = run_model(13,x0,40,'ng/min','SC',false);
[S40_noAGTfb,t40_noAGTfb] = run_model(13,x0,40,'ng/min','SC',false,varargin_noAGTfb);
[S40_noACEfb,t40_noACEfb] = run_model(13,x0,40,'ng/min','SC',false,varargin_noACEfb);

% Indices of key variables
AGT_ind = 1;
PRA_ind = 8;
AngI_circ_ind = 2;
AngI_T_ind = 34;

fold_AGT_sim = S40(AGT_ind,:)/S40(AGT_ind,1);
fold_AGT_sim_noAGTfb = S40_noAGTfb(AGT_ind,:)/S40(AGT_ind,1);

fold_PRA_sim = S40(PRA_ind,:)/S40(PRA_ind,1);
fold_PRA_sim_noAGTfb = S40_noAGTfb(PRA_ind,:)/S40(PRA_ind,1);

fold_AngI_circ_sim = S40(AngI_circ_ind,:)/S40(AngI_circ_ind,1);
fold_AngI_T_sim = S40(AngI_T_ind,:)/S40(AngI_T_ind,1);
fold_AngI_circ_sim_noACEfb = S40_noACEfb(AngI_circ_ind,:)/S40(AngI_circ_ind,1);
fold_AngI_T_sim_noACEfb = S40_noACEfb(AngI_T_ind,:)/S40(AngI_T_ind,1);

%% Figure 

c = summer(6);

figure(4)
subplot(2,3,1)
plot(t40,fold_AGT_sim,'linewidth',1.5,'color','k');
hold on
plot(t40_noAGTfb,fold_AGT_sim_noAGTfb,'--','linewidth',1.5,'color','k');
errorbar(13,fold_AGT,fold_AGT_SD,'o','linewidth',1,...
        'markerfacecolor',c(1,:),'color','k','markersize',10,'linewidth',0.8);
hold off
xlabel('Time (days)');
ylabel('[AGT]_{circ} (ratio to control)');
set(gca,'fontsize',14);
xlim([0,13]);
ttla = title('a','fontsize',20,'fontweight','bold');
ttla.Units = 'Normalize'; 
ttla.Position(1) = 0;
ttla.HorizontalAlignment = 'left';  

subplot(2,3,[2,3])
a1=plot(t40,fold_PRA_sim,'linewidth',1.5,'color','k');
hold on
a2=plot(t40_noAGTfb,fold_PRA_sim_noAGTfb,'--','linewidth',1.5,'color','k');
errorbar([3,7,10],fold_PRA(1:end-1),fold_PRA_SD(1:end-1),'d','linewidth',1,...
        'markerfacecolor',c(1,:),'color','k','markersize',10,'linewidth',0.8);
a3=errorbar(13,fold_PRA(end),fold_PRA_SD(end),'o','linewidth',1,...
        'markerfacecolor',c(1,:),'color','k','markersize',10,'linewidth',0.8);    
hold off
xlabel('Time (days)');
ylabel('PRA (ratio to control)');
legend([a1,a2,a3],{'Simulation','Simulation (fb_{circ(AGT)} = 0)','Data'},'location','northeastoutside');
set(gca,'fontsize',14);
set(legend,'fontsize',14);
xlim([0,13]);
ttl = title('b','fontsize',20,'fontweight','bold');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';  

subplot(2,3,[4,5,6])
b1 = plot(t40,fold_AngI_circ_sim,'linewidth',1.5,'color',c(1,:));
hold on
b2 = plot(t40_noACEfb,fold_AngI_circ_sim_noACEfb,'--','linewidth',1.5,'color',c(1,:));
b3 = errorbar([0,3,7,10],fold_AngI_circ(1:end-1),fold_AngI_circ_SD(1:end-1),'d',...
        'linestyle','none','markerfacecolor',c(1,:),'color','k','markersize',10,'linewidth',0.8);
errorbar(13,fold_AngI_circ(end),fold_AngI_circ_SD(end),'o','linestyle','none',...
    'markerfacecolor',c(1,:),'color','k','markersize',10,'linewidth',0.8);
b4 = plot(t40,fold_AngI_T_sim,'linewidth',1.5,'color',c(4,:));
b5 = plot(t40_noACEfb,fold_AngI_T_sim_noACEfb,'--','linewidth',1.5,'color',c(4,:));
b6 = errorbar([0,3,7,10],fold_AngI_T(1:end-1),fold_AngI_T_SD(1:end-1),'d',...
    'linestyle','none','markerfacecolor',c(4,:),'color','k','markersize',10,'linewidth',0.8);
errorbar(13,fold_AngI_T(end),fold_AngI_T_SD(end),'d','linestyle','none',...
    'markerfacecolor',c(4,:),'color','k','markersize',10,'linewidth',0.8);
hold off
legend([b1,b2,b3,b4,b5,b6],{'[Ang I]_{circ} (Simulation)','[Ang I]_{circ} (Simulation; fb_{circ(ACE)} = 0)','Ang I_{circ} (Data)',...
       '[Ang I]_T (Simulation)','[Ang I]_{T} (Simulation; fb_{circ(ACE)} = 0)','[Ang I]_T (Data)'},'location','northeastoutside');
set(legend,'fontsize',14);
xlabel('Time (days)');
ylabel('[Ang I] (ratio to control)');
set(gca,'fontsize',14);
xlim([0,13]);
ttl = title('c','fontsize',20,'fontweight','bold');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  


