function get_Losartan_fig2()
%% Data
t_vec = [0,0.09,0.14,0.22,0.47,0.76,1,2,4,6,8,10,24]/24; % h 0.76,1,
los_vec = [0,0.69,1.07,1.21,1.26,1.01,0.87,0.78,0.60,0.50,0.39,0.30,0]*10^9/422.91/1000;
los_vec_SD = los_vec - [0,0.65,0.94,1.11,1.14,0.91,0.77,0.68,0.51,0.44,0.31,0.24,0]*10^9/422.91/1000;
t_exp_vec = [0,0.56,0.76,1,2,4,6,8,10,24]/24;
exp_vec = [0,0.13,0.19,0.29,0.48,0.64,0.84,0.97,0.68,0.09]*10^9/436.89/1000;
exp_vec_SD = exp_vec - [0,0.1,0.16,0.2,0.42,0.57,0.8,0.91,0.61,0.07]*10^9/436.89/1000;

%% Simulations
dose = 10*284*10^6/(422.91);
t_res = 1440;

% Run model to get initial Los_GI = dose, compute new IC
SS = load('model_SS.mat').SSdata_los_separate;
num_vars = length(SS); x_p0 = zeros(num_vars,1);
varargin_input1 = {'losartan',dose,'oral; single', true, 12*60/t_res,'all'};
[S,t,YP] = run_sim_renal_RAS(10/24,SS,x_p0,0,'ng/min','none',true,varargin_input1);
SS_new = S(:, end);
x_p1 = YP(:,end);

% Let dose be reabsorbed over 1 day
varargin_input2 = {'losartan',0,'oral; single',false, 10*60/t_res};
[S,t,YP] = run_sim_renal_RAS(1,SS_new,x_p1,0,'ng/min','none',true,varargin_input2);

%% Plot

c = summer(4);
figure(2)
a=plot(t*24, S(76,:),'color','k','linewidth',1.5);
hold on
errorbar(t_vec*24,los_vec,los_vec_SD,'d','color','k','linewidth',0.8,'markerfacecolor','k','markersize',10);
b=plot(t*24, S(78,:),'color',c(2,:),'linewidth',1.5);
errorbar(t_exp_vec*24,exp_vec,exp_vec_SD,'d','color','k','linewidth',0.8,'markerfacecolor',c(2,:),'markersize',10);
hold off
ylabel('Systemic concentration (pmol/mL)')
xlabel('Time (h)')
set(gca,'fontsize',12)
legend([a,b],{'Losartan','EXP3174'})
set(legend,'fontsize',12)
xlim([0,24])

end