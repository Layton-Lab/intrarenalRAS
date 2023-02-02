% Default parameters
    % delta = 0.1;
    % days = 13;
    % D = 40;
    % units = 'ng/min';
    % inf_type = 'SC';
    % D_Los = 30;

function [test, test_los] = get_Losartan_Fig7_and_8(delta,days,D,units,inf_type,D_Los)

separate = true;
t_res = 1440;

% Losartan parameters
los_pars = load('losartan_pars.mat').pars;
los_pars(8) = 0.5; % Delta t
dose = D_Los*284*10^6/(422.91);
   
% Load initial condition
x0 = load('model_SS.mat').SSdata_los_separate;
SS_inf = x0;
num_vars = length(SS_inf);
x_p0_inf = zeros(num_vars,1);
t_end_inf = 0; % to track the time at end of each day
t_tot_inf = []; % to store full time vector
S_tot_inf = []; % to store time series of individual variables
vardata0 = []; % to store time series of combined variables 
avg = []; % to store daily averages of individual variables
avgvars0 = []; % to store daily averages of combined variables


% Get baseline time series (with optimized parameters)
for j = 1:days
    varargin_input1 = {'losartan',dose,'oral; drinking water',false,los_pars(8),'t0',t_end_inf*t_res,'los_pars',los_pars(1:7)};
    [S_inf,t_inf,YP_inf] = run_model(1,SS_inf,x_p0_inf,D,units,inf_type,separate,varargin_input1);
    
    % Update initial conditions
    SS_inf = S_inf(:,end);
    x_p0_inf = YP_inf(:,end);
    
    % Update all vectors
    t_tot_inf = [t_tot_inf, t_inf + t_end_inf]; 
    t_end_inf = t_tot_inf(end);
    S_tot_inf = [S_tot_inf S_inf]; 
    avg = [avg, mean(S_inf,2)];
    vars = compute_regional_concentrations(S_inf);
    vardata0 = [vardata0 vars];
    avgvars0 = [avgvars0, mean(vars,2)];
end

% Store baseline absolute and average variable values at end of experiment
SS_0 = S_tot_inf(:,end);
SS_avg_0 = avg(:,end);

% Initialize variables to store sensitivity analysis results
p_index = 8; % Number of parameters
SSdata = zeros(num_vars,p_index); % to store absolute values of individual vars at end of experiment
avgdata = zeros(num_vars,p_index); % to store average values of individual vars on final day
sens_SS = zeros(num_vars,p_index); % to store fold-change in absolute values of individual vars at end of experiment
sens_avg_SS = zeros(num_vars,p_index); % to store fold-change in average values of individual vars on final day
vardata = zeros(12,p_index); % to store absolute values of combined vars at end of experiment
avgvardata = zeros(12,p_index); % to store average values of combined vars on final day

for i = 1:p_index

     % Perturb parameter i by delta
     pars_i = los_pars;
     delta_pi = delta*los_pars;
     pars_i(i) = pars_i(i) + delta_pi(i);
     
     % Re-initialize variables
     SS_inf = x0;
     x_p0_inf = zeros(num_vars,1);
     t_end_inf = 0;
     t_tot_inf = [];
     S_tot_inf = [];
     avg = [];
     new_vars = [];
     avg_vars = [];

     % Simulate infusion + Losartan for desired number of days
     for j = 1:days

        varargin_input1 = {'losartan',dose,'oral; drinking water',false,pars_i(8),'t0',t_end_inf*1440,'los_pars',pars_i(1:7)};
        [S_inf,t_inf,YP_inf] = run_model(1,SS_inf,x_p0_inf,D,units,inf_type,separate,varargin_input1);
    
        % Update initial conditions
        x_p0_inf = YP_inf(:,end);
        SS_inf = S_inf(:,end);

        % Update all vectors
        t_tot_inf = [t_tot_inf, t_inf + t_end_inf];
        t_end_inf = t_tot_inf(end);
        S_tot_inf = [S_tot_inf S_inf(:,:)];
        avg = [avg, mean(S_inf(:,:),2)];
        vars = compute_regional_concentrations(S_inf);
        new_vars = [new_vars vars];
        avg_vars = [avg_vars, mean(vars,2)];

        j
     end

     % Store results for each parameter perturbation
     sim = S_tot_inf;
     SSdata(:,i) = sim(:,end);
     avgdata(:,i) = avg(:,end);
     sens_SS(:,i) = (SSdata(:,i) - SS_0)./SS_0;
     sens_avg_SS(:,i) = (avgdata(:,i) - SS_avg_0)./SS_avg_0;
     vardata(:,i) = new_vars(:,end);
     avgvardata(:,i) = avg_vars(:,end);

     i
end

%% Get plotting results for the effect on RAS variables

% Compare average variable values on the final day of the experiment
res = [avgvardata;avgdata(1:6,:);avgdata(8,:);avgdata(58,:);avgdata(59,:);avgdata(60,:)]; % perturbed
res0 = [avgvars0(:,end);SS_avg_0(1:6);SS_avg_0(8);SS_avg_0(58);SS_avg_0(59);SS_avg_0(60)]; % control

sens = (res - res0)./res0;

% Any variable change < 5% is ignored to improve plot readability
delta_SS = 100*sens;
test = delta_SS;
[r,c] = find(abs(delta_SS) <= 5);
for i = 1:length(r)
    test(r(i),c(i)) = NaN;
end

% Extract desired variables for plotting
test_circ = [test(13,:);test(19,:);test(14,:);test(15,:);test(end-2,:);test(end-1,:);test(end,:)];
test_intra = test(1:3,:);
test_memb = test(4:6,:);
test_ext = test(7:9,:);
test_WK = test(10:12,:);

%% Plot sensitivity results of RAS variables

params = {'k_a^{Los}','k_{cyt}', 'k_{sp}', 'k_{ps}','k_{elim}^{Los}','k_{elim}^{EXP3174}', 'B_{AT1R}^{+}', '\Deltat'};

circRowNames = {'[AGT]_{circ}', 'PRA', '[Ang I]_{circ}','[Ang II-p]_{circ}',...
            '[Ang II-i]_{circ}', '[Ang II]_{circ}','[AT1R-bound Ang II]_{circ}'};

intraRowNames = {'Endogenous', 'Exogenous' 'Total'}; %'[Ang II-p]^{Cell}_{T}', '[Ang II-i]^{Cell}_{T}','[Ang II]^{Cell}_{T}'};

membRowNames = {'Endogenous', 'Exogenous' 'Total'};
%{'[AT1R-bound Ang II-p]^{Memb}_{T}', '[Ang II-i]^{Memb}_{T}','[Ang II]^{Memb}_{T}'};

extRowNames = {'Endogenous', 'Exogenous' 'Total'};
%{'[Ang II-p]^{Ext}_{T}', '[Ang II-i]^{Ext}_{T}','[Ang II]^{Ext}_{T}'};

WKRowNames = {'Endogenous', 'Exogenous' 'Total'};
%{'[Ang II-p]_{T}', '[Ang II-i]_{T}','[Ang II]_{T}'};

figure(8)
subplot(5,1,1)
heatmap(params,circRowNames,test_circ,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<5%',...
    'ColorLimits',[min(test,[],'all'),max(test,[],'all')],...
    'ColorbarVisible','on');
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(5,1,2)
heatmap(params,membRowNames,test_memb,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<5%',...%,...
       'ColorLimits',[min(test,[],'all'),max(test,[],'all')],...
    'ColorbarVisible','off');
ylabel({'Membrane- and', 'AT1R-bound [Ang II]'})
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(5,1,3)
heatmap(params,intraRowNames,test_intra,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<5%',...
        'ColorLimits',[min(test,[],'all'),max(test,[],'all')],'ColorbarVisible','off',...
        'ColorbarVisible','off');
ylabel({'Intracellular', '[Ang II]'})
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(5,1,4)
heatmap(params,extRowNames,test_ext,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<5%',...
       'ColorLimits',[min(test,[],'all'),max(test,[],'all')],...
    'ColorbarVisible','off');
ylabel({'Extracellular', '[Ang II]'})
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(5,1,5)
heatmap(params,WKRowNames,test_WK,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<5%',...,...
        'ColorLimits',[min(test,[],'all'),max(test,[],'all')],...
    'ColorbarVisible','off');  
ylabel({'Whole', 'kidney [Ang II]'})
xlabel('Parameters');
set(gca,'fontsize',10);

%% Get plotting results for the effect on the Losartan-EXP3174 variables

delta_SS = 100*sens_avg_SS; % sens_SS
test_los = delta_SS;
[r,c] = find(abs(delta_SS) <= 1);
for i = 1:length(r)
    test_los(r(i),c(i)) = NaN;
end

los_i = 76;

test_sys = test_los(los_i:los_i+3,:);
test_renal = test_los(los_i+24:los_i+25,:);
test_per = test_los(los_i+6:los_i+7,:);
test_GI = test_los(los_i+4,:);

%% Plot sensitivity results of Losartan variables

sysRowNames = {'[Los]_{circ}', '[AT1R-bound Los]_{circ}', ...
               '[EXP3174]_{circ}', '[AT1R-bound EXP3174]_{circ}' };

renalRowNames = {'[Los]_{T}', '[EXP3174]_T'};

perRowNames = {'[Los]_{peri}', '[EXP3174]_{peri}'};

GIRowNames = {'[Los]_{GI}'};

figure(7)
subplot(4,1,1)
heatmap(params,sysRowNames,test_sys,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<5%',...
    'ColorLimits',[min(test_los(los_i:los_i+25,:),[],'all'),max(test_los(los_i:los_i+25,:),[],'all')],...
    'ColorbarVisible','on');
%ylabel({'Systemic', 'Compartment'})
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(4,1,2)
heatmap(params,renalRowNames,test_renal,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...%,...
       'ColorLimits',[min(test_los(los_i:los_i+25,:),[],'all'),max(test_los(los_i:los_i+25,:),[],'all')],...
    'ColorbarVisible','off');
%ylabel({'Renal', 'Compartment'})
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(4,1,3)
heatmap(params,perRowNames,test_per,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test_los(los_i:los_i+25,:),[],'all'),max(test_los(los_i:los_i+25,:),[],'all')],'ColorbarVisible','off',...
        'ColorbarVisible','off');
%ylabel({'Peripheral', 'Compartment'})
set(gca,'fontsize',10);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));


subplot(4,1,4)
heatmap(params,GIRowNames,test_GI,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
       'ColorLimits',[min(test_los(los_i:los_i+25,:),[],'all'),max(test_los(los_i:los_i+25,:),[],'all')],...
    'ColorbarVisible','off');
%ylabel({'Gastrointestinal', 'Compartment'})
set(gca,'fontsize',10);
xlabel('Parameters');


end
