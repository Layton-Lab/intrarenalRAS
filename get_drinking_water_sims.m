% Default parameters
    % D_Los = 30 (mg/kg)
    % separate = true

function [t_tot, S_tot, avgs, t_tot_inf, S_tot_inf, avgs_inf] = get_drinking_water_sims(D_Los, separate)

deltat = 0.5;

dose = D_Los*284*10^6/(422.91);
varargin_input0 = {'losartan',dose,'oral; drinking water',false,deltat};

if separate
    SS = load('model_SS.mat').SSdata_los_separate;
    los_i = 76; % First index corresponding to Losartan variables
else
    SS = load('model_SS.mat').SSdata_los; 
    los_i = 42;
end

% Initialize all variables needed to solve the DE
SS_inf = SS;
num_vars = length(SS);
x_p0 = zeros(num_vars,1); 
x_p0_inf = zeros(num_vars,1); 
S_tot = []; S_tot_inf = []; 
t_tot = []; t_tot_inf = [];
t_end = 0; t_end_inf = 0;

% Simulate 13 days of Losartan administration in drinking water
for i = 1:13 

    varargin_input1 = {'losartan',dose,'oral; drinking water',false,deltat,'t0',t_end_inf*1440};

    % Control
    [S,t,YP] = run_model(1,SS,x_p0,0,'ng/min','none',separate,varargin_input0);

    x_p0 = YP(:,end);
    SS = S(:,end);
    S_tot = [S_tot S];
    t_tot = [t_tot, t + t_end];
    t_end = t_tot(end);
    

    % Ang II infusion
    [S_inf,t_inf,YP_inf] = run_model(1,SS_inf,x_p0_inf,40,'ng/min','SC',separate,varargin_input1);
    
    x_p0_inf = YP_inf(:,end);
    SS_inf = S_inf(:,end);
    S_tot_inf = [S_tot_inf S_inf];
    t_tot_inf = [t_tot_inf, t_inf + t_end_inf];
    t_end_inf = t_tot_inf(end);

    i
end

avgs = [];
avgs_inf = [];
for i = 0:12
    [~,t_ind0] = min(abs(t_tot-i*ones(size(t_tot))));
    [~,t_ind1] = min(abs(t_tot-(i+1)*ones(size(t_tot))));
    [~,t_ind0_inf] = min(abs(t_tot_inf-i*ones(size(t_tot_inf))));
    [~,t_ind1_inf] = min(abs(t_tot_inf-(i+1)*ones(size(t_tot_inf))));

    dayi = S_tot(:,t_ind0:t_ind1);
    avg_dayi = mean(dayi, 2);
    avgs = [avgs,avg_dayi];
 
    dayi_inf = S_tot_inf(:,t_ind0_inf:t_ind1_inf);
    avg_dayi_inf = mean(dayi_inf, 2); 
    avgs_inf = [avgs_inf,avg_dayi_inf];

end

end