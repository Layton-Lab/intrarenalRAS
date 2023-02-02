% Default parameters
    % days = 13;
    % D = 40;
    % units = 'ng/min';
    % inf_type = 'SC';
    % D_Los = 30;

function get_Losartan_Fig6(days,D,units,inf_type,D_Los, days_Los)

deltat = 0.5;
t_res = 1440;

dose = D_Los*284*10^6/(422.91);

% Run infusion for 2 weeks, without Losartan administration
SS = load('model_SS.mat').SSdata_separate;
x_p0 = zeros(length(SS),1);
[S40,t40,YP40] = run_sim_renal_RAS(days,SS,x_p0,D,units,inf_type,true);

% Take end both of previous simulation as the new initial condtion
SS_inf = [S40(:,end);zeros(26,1)];
x_p0_inf = [YP40(:,end);zeros(26,1)]; 
S_tot_inf = [S40; zeros(26,length(t40))]; t_tot_inf = t40; 
t_end_inf = t_tot_inf(end);

% Simulation Losartan administration in drinking water for 1 week
for i = 0:(days_Los-1)

    varargin_input1 = {'losartan',dose,'oral; drinking water',false,deltat,'t0',t_end_inf*t_res};

    % Ang II infusion, mechanism (i)
    [S_inf,t_inf,YP_inf] = run_model(1,SS_inf,x_p0_inf,40,'ng/min','SC',true,varargin_input1);
    
    x_p0_inf = YP_inf(:,end);
    SS_inf = S_inf(:,end);
    S_tot_inf = [S_tot_inf S_inf(:,:)];
    t_tot_inf = [t_tot_inf, t_inf + t_end_inf];
    t_end_inf = t_tot_inf(end);
    
    i
end

% Compute exogenous vs. endogenous distribution
t = t_tot_inf;
S = S_tot_inf;
[comp_intra_tot_los, comp_memb_tot_los, comp_isf_tot_los, comp_tot_los] = get_endo_exo_distributions(S);

%% Data
AngII_T_con = [164,132,133,138,161];
AngII_T = [164,148,165,209,350];
AngII_T_SD = [23,15,8,10,62];
fold_AngII_T = AngII_T./AngII_T_con;
fold_AngII_T_SD = AngII_T_SD./AngII_T_con;

%% Plot

c = summer(5);

figure(6)

subplot(2,3,[1,4])
errorbar([0,3,7,10,13],fold_AngII_T*S_tot_inf(75,1),fold_AngII_T_SD*S_tot_inf(75,1),'d','color','k',...
    'linewidth',1,'markerfacecolor','k','markersize',8);
hold on
plot(t_tot_inf,S_tot_inf(75,:),'color','k','LineWidth',1);
b=plot(t_tot_inf,S_tot_inf(35,:),'color',c(3,:),'Linewidth',1);
d=plot(t_tot_inf,S_tot_inf(58,:),'color',c(1,:),'Linewidth',1);
xline(13,':','color','k','linewidth',1.5)
hold off
set(gca,'fontsize',12)
xlabel('Time (days)')
ylabel('Whole kidney [Ang II] (fmol/g kidney)')
xlim([0,20])
ttla = title('A', 'fontsize',14);
ttla.Units = 'Normalize'; 
ttla.Position(1) = 0;
ttla.HorizontalAlignment = 'left'; 
legend([b,d],{'Endogenous','Exogenous'},'fontsize',11)
hold off

subplot(2,3,2)
area(t,comp_tot_los(:,1)+comp_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
set(gca,'Fontsize',12);
ylabel({'Fraction of whole', 'kidney Ang II'});
xlabel('Time (days)');
ylim([0,1.1]);
xlim([0,20]);
ttlb = title('B', 'fontsize',14);
ttlb.Units = 'Normalize'; 
ttlb.Position(1) = 0;
ttlb.HorizontalAlignment = 'left'; 
hold off

subplot(2,3,3)
area(t,comp_memb_tot_los(:,1)+comp_memb_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_memb_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
set(gca,'Fontsize',12);
ylabel({'Membrane-bound fraction', 'of whole kidney Ang II'});
xlabel('Time (days)');
xlim([0,20]);
ttld = title('C', 'fontsize',14);
ttld.Units = 'Normalize'; 
ttld.Position(1) = 0;
ttld.HorizontalAlignment = 'left'; 
hold off
ylim([0,0.1])


subplot(2,3,5)
area(t,comp_intra_tot_los(:,1)+comp_intra_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_intra_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
set(gca,'Fontsize',12);
ylabel({'Intracellular fraction of', 'whole kidney Ang II'});
xlabel('Time (days)');
ttlc = title('D', 'fontsize',14);
ttlc.Units = 'Normalize'; 
ttlc.Position(1) = 0;
ttlc.HorizontalAlignment = 'left'; 
xlim([0,20])
ylim([0,0.6])

subplot(2,3,6)
area(t,comp_isf_tot_los(:,1) + comp_isf_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_isf_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
set(gca,'Fontsize',12);
ylabel({'Extracellular fraction of', 'whole kidney Ang II'});
xlabel('Time (days)');
ylim([0,1.1]);
xlim([0,20]);
ttld = title('E', 'fontsize',14);
ttld.Units = 'Normalize'; 
ttld.Position(1) = 0;
ttld.HorizontalAlignment = 'left'; 

end

