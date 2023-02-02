function get_Losartan_Fig3_to_5()
%% Ang II infusion simulation (experiment (i))

x0 = load('model_SS.mat').SSdata_separate;
xp0 = zeros(length(x0),1);
[Si,ti] = run_model(13,x0, xp0, 40,'ng/min','SC',true);

% AGT, PRA, AngI, AngII, AngI_T, AngII_T
con = [x0(1); x0(8); x0(2); x0(59); x0(34); x0(75)]; 
AngIIi = [Si(1,end); Si(8,end); Si(2,end); Si(59,end); Si(34,end); Si(75,end)];

% Get daily averages 
avgs_40 = [];
for j = 0:12
    [~,t_ind0_40] = min(abs(ti-j*ones(size(ti))));
    [~,t_ind1_40] = min(abs(ti-(j+1)*ones(size(ti))));
 
    dayj = Si(:,t_ind0_40:t_ind1_40);
    avg_dayj = mean(dayj, 2); 
    avgs_40 = [avgs_40,avg_dayj];
end

%% Losartan simulations (experiment (ii) and (iii))
[~, S_tot, avgs, t_tot_inf, S_tot_inf, avgs_inf] = get_drinking_water_sims(30, true);

AGT = avgs(1,end); AGT_inf = avgs_inf(1,end); 
PRA = avgs(8, end); PRA_inf = avgs_inf(8, end); 
AngI = avgs(2, end); AngI_inf = avgs_inf(2, end); 
AngII = avgs(59, end); AngII_inf = avgs_inf(59, end); 
AngI_T = avgs(34, end); AngI_T_inf = avgs_inf(34, end);
AngII_T = avgs(75, end); AngII_T_inf = avgs_inf(75, end);

los_con = [AGT; PRA; AngI; AngII; AngI_T; AngII_T];
los_AngIIi = [AGT_inf; PRA_inf; AngI_inf; AngII_inf; AngI_T_inf; AngII_T_inf];

%% Data
AngIID = [1.58; 0.05/4.17; 0.165; 2.936; 104/122; 350/161];
AngIIS = [1.75-1.58; 0.02/4.17; 0.215-0.165; 3.383 - 2.936; (120-104)/122; 62/161];

los_conD = [0.65; 11.52; 7.87; 5.83; 1.09; 0.75];
los_conS = [0.14; 0.54; 1.25; 0.48; 0.25; 0.06];
los_AngIID = [295.9/498.4; 50.44/4.17; 1566.6/213.1; 289.6/50.4; 258.1/121.2; 183.7/165.2];
los_AngIIS = [(344.1-295.9)/498.4; 8.33/4.17; (1825.9-1566.6)/213.1; (341.7-289.6)/50.4; (298.3-258.1)/121.2; (207.8 - 183.7)/165.2];

%% Losartan model validation (Figure 3)
names3 = {'[AGT]_{circ}','PRA','[AngI]_{circ}','[AngII]_{circ}','[AngI]_T','[AngII]_T'};

b1 = [0.75-0.125, 3.75-0.125, 6.75-0.125, 9.75-0.125, 12.75-0.125, 15.75-0.125];

b2 = [1.25-0.125, 4.25-0.125, 7.25-0.125, 10.25-0.125, 13.25-0.125, 16.25-0.125];

b3 = [1.75+0.125, 4.75+0.125, 7.75+0.125, 10.75+0.125, 13.75+0.125, 16.75+0.125];

b4 = [2.25+0.125, 5.25+0.125, 8.25+0.125, 11.25+0.125, 14.25+0.125, 17.25+0.125];

c = summer(16);

figure(3)
bar(b1, con,'BarWidth', 0.15, 'facecolor', c(1,:), 'facealpha', 0.4)
hold on
bar(b3, los_con,'BarWidth', 0.15, 'facecolor', c(11,:), 'facealpha', 0.4)
b=bar(b1, con,'barwidth', 0.1,'facecolor', c(1,:));
l=bar(b3, los_conD.*con,'barwidth', 0.1, 'facecolor', c(11,:));
errorbar(b3, los_conD.*con, los_conS.*con, 'linestyle','none','color','k','linewidth',0.8);

bar(b2, AngIIi,'barwidth', 0.15, 'facecolor', c(6,:), 'facealpha', 0.4)
bar(b4, los_AngIIi,'barwidth', 0.15, 'facecolor', c(16,:), 'facealpha', 0.4)
g=bar(b2, AngIID.*con,'barwidth', 0.1,'facecolor', c(6,:));
h=bar(b4, los_AngIID.*con,'barwidth', 0.1, 'facecolor', c(16,:));
errorbar(b2, AngIID.*con, AngIIS.*con, 'linestyle','none','color','k','linewidth',0.8);
errorbar(b4, los_AngIID.*con, los_AngIIS.*con, 'linestyle','none','color','k','linewidth',0.8);

hold off
legend([b,g,l,h],{'Control', '(i) Ang II', '(ii) Losartan','(iii) Ang II + Losartan'}, 'fontsize', 12)
set(gca, 'xtick', [1.5, 4.5, 7.5, 10.5, 13.5 16.5],'xticklabel', names3, 'yscale','log', 'fontsize', 12)
xlim([0,18])
ylim([5,5*10^6])
ylabel({'Log concentrtion (fmol/mL)', 'or activity (fmol/mL per min)'})

%%  Comparing plasma and intrarenal Ang II accumulation (Figure 4)

c = summer(5);

figure(4)
subplot(1,2,1)
box on
hold on
e1 = bar(1,[S_tot(3,1),0],'stacked', 'Facecolor','flat');
e1(1).CData = c(3,:);
e1(2).CData = c(1,:);

d1 = bar(2,[avgs_40(3,end),avgs_40(42,end)],'stacked', 'Facecolor','flat');
d1(1).CData = c(3,:);
d1(2).CData = c(1,:);

c4 = bar(3,[avgs(3,end),avgs(42,end)],'stacked', 'Facecolor','flat');
c4(1).CData = c(3,:);
c4(2).CData = c(1,:);

c5 = bar(4, [avgs_inf(3,end),avgs_inf(42,end)],'stacked','Facecolor','flat');
c5(1).CData = c(3,:);
c5(2).CData = c(1,:);
set(gca,'xtick',1:4,'xticklabel',...
    {'Control','(i) Ang II', '(ii) Losartan','(iii) Ang II + Losartan'}, 'fontsize', 12)
ylabel('[Ang II]_{circ} (fmol/mL)')
ttla = title('A', 'fontsize',14);
ttla.Units = 'Normalize'; 
ttla.Position(1) = 0;
ttla.HorizontalAlignment = 'left'; 
hold off

subplot(1,2,2)
e2 = bar(1,[S_tot(35,1),0],'stacked', 'Facecolor','flat');
e2(1).CData = c(3,:);
e2(2).CData = c(1,:);
hold on

d2 = bar(2,[avgs_40(35,end),avgs_40(58,end)],'stacked', 'Facecolor','flat');
d2(1).CData = c(3,:);
d2(2).CData = c(1,:);
hold on
c1 = bar(3,[avgs(35,end),avgs(58,end)],'stacked','FaceColor','flat');
c1(1).CData = c(3,:);
c1(2).CData = c(1,:);
c2 = bar(4, [avgs_inf(35,end),avgs_inf(58,end)],'stacked', 'FaceColor','flat');
c2(1).CData = c(3,:);
c2(2).CData = c(1,:);
hold off
set(gca,'xtick',1:4,'xticklabel',...
    {'Control','(i) Ang II', '(ii) Losartan','(iii) Ang II + Losartan'}, 'fontsize', 12)
legend('Endogenous', 'Exogenous');
set(legend, 'fontsize',12);
ylabel('[Ang II]_T (fmol/mL)')
ttlb = title('B', 'fontsize',14);
ttlb.Units = 'Normalize'; 
ttlb.Position(1) = 0;
ttlb.HorizontalAlignment = 'left'; 

%% Intrarenal distribution of endogenous and exogenous Ang II (Figure 5)
[comp_intra_tot, comp_memb_tot, comp_isf_tot, comp_tot] = get_endo_exo_distributions(Si);
[comp_intra_tot_los, comp_memb_tot_los, comp_isf_tot_los, comp_tot_los] = get_endo_exo_distributions(S_tot_inf);

t = t_tot_inf;

c = summer(5);

figure(5)

subplot(2,4,1) % Total renal Ang II
area(ti,comp_tot(:,1)+comp_tot(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(ti,comp_tot(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttla=title('A.i. Ang II','fontsize',10);
set(gca,'Fontsize',12,'xtick',[0,5,10]);
ylabel('Fraction of whole kidney Ang II');
xlabel('Time (days)')
ylim([0,1.1]);
xlim([0,13]);
ttla.Units = 'Normalize'; 
ttla.Position(1) = 0;
ttla.HorizontalAlignment = 'left';  

subplot(2,4,2) % Total renal Ang II with los
area(t,comp_tot_los(:,1)+comp_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttlb=title('A.iii. Ang II + Losartan','fontsize',10);
set(gca,'Fontsize',10, 'yticklabel',[],'xtick',[0,5,10]);
xlabel('Time (days)')
ylim([0,1.1]);
xlim([0,13]);
ttlb.Units = 'Normalize'; 
ttlb.Position(1) = 0;
ttlb.HorizontalAlignment = 'left'; 

subplot(2,4,3) % Membrane bound Ang II
area(ti,comp_memb_tot(:,1)+comp_memb_tot(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(ti,comp_memb_tot(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttle= title('B.i. Ang II','fontsize',10);
set(gca,'Fontsize',12);
ylabel({'Membrane-bound fraction','of whole kidney Ang II'});
xlabel('Time (days)');
ylim([0,0.1]);
xlim([0,13]);
ttle.Units = 'Normalize'; 
ttle.Position(1) = 0;
ttle.HorizontalAlignment = 'left';

subplot(2,4,4) % Membrane bound Ang II with los
area(t,comp_memb_tot_los(:,1)+comp_memb_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_memb_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttlf = title('B.iii. Ang II + Losartan','fontsize',10);
set(gca,'Fontsize',10, 'yticklabel',[]);
xlabel('Time (days)');
ylim([0,0.1]);
xlim([0,13]);
ttlf.Units = 'Normalize'; 
ttlf.Position(1) = 0;
ttlf.HorizontalAlignment = 'left';

axes('Position',[0.15,0.7,0.15,0.2]);
hold on
box on
area(t*24,comp_memb_tot_los(:,1)+comp_memb_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
area(t*24,comp_memb_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
xlim([0,5])
set(gca, 'YTickLabel',[],'xtick',[0,1,2,3,4,5])
ylim([0,0.1])
xlabel('Time (hours)')

subplot(2,4,5) % Intracellular Ang II
area(ti,comp_intra_tot(:,1)+comp_intra_tot(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(ti,comp_intra_tot(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttlc= title('C.i. Ang II','fontsize',10);
set(gca,'Fontsize',10);
ylabel({'Intracellular fraction','of whole kidney Ang II'});
xlabel('Time (days)');
ylim([0,0.6]);
xlim([0,13]);
ttlc.Units = 'Normalize'; 
ttlc.Position(1) = 0;
ttlc.HorizontalAlignment = 'left';

subplot(2,4,6) % Intracellular Ang II with los
area(t,comp_intra_tot_los(:,1)+comp_intra_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_intra_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttld = title('C.iii. Ang II + Losartan','fontsize',10);
set(gca,'Fontsize',10, 'yticklabel',[]);
xlabel('Time (days)');
ylim([0,0.6]);
xlim([0,13]);
ttld.Units = 'Normalize'; 
ttld.Position(1) = 0;
ttld.HorizontalAlignment = 'left';

axes('Position',[0.15,0.7,0.15,0.2]);
hold on
box on
area(t*24,comp_intra_tot_los(:,1)+comp_intra_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
area(t*24,comp_intra_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
xlim([0,5])
set(gca, 'YTickLabel',[],'xtick',[0,1,2,3,4,5])
ylim([0,0.6])
xlabel('Time (hours)')

subplot(2,4,7) % Total renal Ang II
area(ti,comp_isf_tot(:,1)+comp_isf_tot(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(ti,comp_isf_tot(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
ttlg= title('D.i. Ang II','fontsize',10);
set(gca,'Fontsize',10);
ylabel({'Extracellular fraction','of whole kidney Ang II'});
xlabel('Time (days)');
ylim([0,1.1]);
xlim([0,13]);
ttlg.Units = 'Normalize'; 
ttlg.Position(1) = 0;
ttlg.HorizontalAlignment = 'left';

subplot(2,4,8) % Total renal Ang II with los
area(t,comp_isf_tot_los(:,1)+comp_isf_tot_los(:,2),'facecolor',c(1,:));%,'linewidth',1.5);
hold on
area(t,comp_isf_tot_los(:,1),'facecolor',c(3,:));%,'linewidth',1.5);
hold off
legend('Exogenous', 'Endogenous')
set(legend,'fontsize',10)
ttlh = title('D.iii. Ang II + Losartan','fontsize',10);
set(gca,'Fontsize',10, 'yticklabel',[]);
xlabel('Time (days)');
ylim([0,1.1]);
xlim([0,13]);
ttlh.Units = 'Normalize'; 
ttlh.Position(1) = 0;
ttlh.HorizontalAlignment = 'left';

end