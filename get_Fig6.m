% Plots losartan distribution figure
function get_Fig6(dose,unit,type,days)

separate = false;
x0 = [load('model_SS.mat').SSdata];
[S_tot,t_tot] = run_model(days,x0,dose,unit,type,separate);

[comp_intra_tot, comp_memb_tot, comp_isf_tot, comp_tot] = get_distributions(S_tot,separate);

c = summer(8);

figure(6)

subplot(2,3,[1,4])
area((t_tot),comp_tot(:,1)+comp_tot(:,2)+comp_tot(:,3)+comp_tot(:,4),'facecolor',c(1,:),'linewidth',1.5);
hold on
area((t_tot),comp_tot(:,1)+comp_tot(:,2)+comp_tot(:,3),'facecolor',c(3,:));
area((t_tot),comp_tot(:,1)+comp_tot(:,2),'facecolor',c(5,:));
area((t_tot),comp_tot(:,1),'facecolor',c(7,:));
hold off
ttl = title('A','fontsize',20);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
set(gca,'Fontsize',16);
ylabel('Fraction of whole kidney Ang II');
xlabel('Time (days)');
ylim([0,1.1]);
xlim([0,13]);

subplot(2,3,2)
area((t_tot),comp_isf_tot(:,1)+comp_isf_tot(:,2)+comp_isf_tot(:,3)+comp_isf_tot(:,4),'facecolor',c(1,:),'linewidth',1.5);
hold on
area((t_tot),comp_isf_tot(:,1)+comp_isf_tot(:,2)+comp_isf_tot(:,3),'facecolor',c(3,:));
area((t_tot),comp_isf_tot(:,1)+comp_isf_tot(:,2),'facecolor',c(5,:));
area((t_tot),comp_isf_tot(:,1),'facecolor',c(7,:));
hold off
ttl = title('B','fontsize',20);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';  
set(gca,'Fontsize',16);
ylabel({'Extracellular fraction of', 'whole kidney Ang II'});
xlabel('Time (days)');
xlim([0,13]);
ylim([0,0.7]);

subplot(2,3,3)
area((t_tot),comp_memb_tot(:,1)+comp_memb_tot(:,2)+comp_memb_tot(:,3)+comp_memb_tot(:,4),'facecolor',c(1,:),'linewidth',1.5);
hold on
area((t_tot),comp_memb_tot(:,1)+comp_memb_tot(:,2)+comp_memb_tot(:,3),'facecolor',c(3,:));
area((t_tot),comp_memb_tot(:,1)+comp_memb_tot(:,2),'facecolor',c(5,:));
area((t_tot),comp_memb_tot(:,1),'facecolor',c(7,:));
hold off
ttl = title('C','fontsize',20);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
set(gca,'Fontsize',16);
ylabel({'Membrane-bound fraction', 'of whole kidney Ang II'});
xlabel('Time (days)');
xlim([0,13]);

subplot(2,3,[5,6])
area((t_tot),comp_intra_tot(:,1)+comp_intra_tot(:,2)+comp_intra_tot(:,3)+comp_intra_tot(:,4),'facecolor',c(1,:),'linewidth',1.5);
hold on
area((t_tot),comp_intra_tot(:,1)+comp_intra_tot(:,2)+comp_intra_tot(:,3),'facecolor',c(3,:));
area((t_tot),comp_intra_tot(:,1)+comp_intra_tot(:,2),'facecolor',c(5,:));
area((t_tot),comp_intra_tot(:,1),'facecolor',c(7,:));
legend({'Renal vasculature compartment','Tubular compartment',...
        'Peritubular compartment','Glomerular compartment'},...
        'location','eastoutside');
hold off
ttl = title('D','fontsize',20);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
set(gca,'Fontsize',16);
ylabel({'Intracellular fraction of', 'whole kidney Ang II'});
xlabel('Time (days)');
xlim([0,13]);
set(legend,'fontsize',16);

end