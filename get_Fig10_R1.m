function [SSdata,test] = get_Fig10_R1(delta,days,D,units,inf_type)
% Resolution (number of minutes/hour)
t_res = 1440; 
    
% Number of points for plotting resolution
N = ((days)*t_res);
    
% Initial time (min); Final time (min); Time vector
t0 = 0*t_res; tend = days*t_res; tspan = linspace(t0,tend,N);

% Load original parameters
pars = get_pars;

% Don't infusion Ang II
infusion_SS.dose = 0; infusion_SS.type = 'none';

% Don't separate exogenous from endogenous Ang II
infusion_SS.separate = false;

% Turn off all feedback
varargin_input_SS = {'all','renin'}; 
  
% Scale infusion dose based on units
if strcmp(units,'ng/kg/min')
    infusion.dose = D*955.84*0.284;
elseif strcmp(units,'ng/min')
    infusion.dose = D*955.84;
end
    
infusion.type = inf_type;
infusion.separate = false;
   
% Load initial condition
x0 = [load('model_SS.mat').SSdata];
num_vars = length(x0); x_p0 = zeros(num_vars,1);

% Set solver options
options = odeset('RelTol',1e-6,'AbsTol',1e-4);

% Simulate infusion with original parameters
sol_orig_pars = ode15i(@(t,x,x_p) ...
                model(t,x,x_p,pars,infusion), ...
                            tspan, x0, x_p0,options);  

sim_0 = sol_orig_pars.y;
SS_0 = sim_0(:,end);

p_index = 15 + 8;
SSdata = zeros(num_vars,p_index);
sens_SS = zeros(num_vars,p_index);

for i = 1:p_index   
     pars_i = pars;
     delta_pi = delta*pars(i+19);
     pars_i(i+19) = pars(i+19) + delta_pi;
     
     % Solve for new initial condition with parameter set pars_i
     sol_0 = ode15i(@(t,x,x_p) ...
                model(t,x,x_p,pars_i,infusion_SS,varargin_input_SS{:}), ...
                            tspan, x0, x_p0,options);  
     S_0 = sol_0.y;
     
     x0_new = S_0(:,end);
     
     varargin_input = {'new_SS',x0_new(11),x0_new(18),x0_new(25),x0_new(4)}; 
     sol = ode15i(@(t,x,x_p) ...
                model(t,x,x_p,pars_i,infusion,varargin_input{:}), ...
                            tspan, x0_new, x_p0,options);  
     sim = sol.y;
     SSdata(:,i) = sim(:,end);
     sens_SS(:,i) = (SSdata(:,i) - SS_0)./SS_0;
end

%% Plotting

params = {'k_{int}','k_{rec}', 'k_{lys}','k_{meg}', 'k_{sec}','k_{diff}', 'AT1R_{Gl}^{tot}','c_{chym}','c_{NEP}',...
          'v_{max}','k_{AngI_{Pt}}','c_{ACE_{Pt}}','k_{AngI_{Tb}}','c_{ACE_{Tb}}',...
          'k_{AGT}','k_a','B_{AT1R}^{-}','K_{circ_{ACE}}','K_{circ_{AGT}}',...
          'K_{Pt}','K_{Tb}','AT1R_{circ_{tot}}','AT1R_{Pv_{tot}}'};
      
% Preferred Order of SS variables (match table in paper)
circRowNames = {'[AGT]_{circ}', '[Ang I]_{circ}','[Ang II]_{circ}',...
            '[AT1R-bound Ang II]_{circ}^{Memb}', '[AT1R]_{circ}^{Memb}'...
            '[Ang(1-7)]_{circ}','PRA'};
        
GlRowNames = {'[Ang I]_{Gl}^{Isf}','[Ang II]_{Gl}^{Isf}','[AT1R-bound Ang II]_{Gl}^{Memb}',...
               '[AT1R-bound Ang II]_{Gl}^{Cell}','[Ang II]_{Gl}^{Cell}',...
               '[AT1R]_{Gl}^{Memb}','[AT1R]_{Gl}^{Cell}'};
           
PtRowNames = {'[Ang I]_{Pt}^{Isf}' ,'[Ang II]_{Pt}^{Isf}','[AT1R-bound Ang II]_{Pt}^{Memb}',...
              '[AT1R-bound Ang II]_{Pt}^{Cell}','[Ang II]_{Pt}^{Cell}',...
              '[AT1R]_{Pt}^{Memb}','[AT1R]_{Pt}^{Cell}'}; 
          
TbRowNames = {'[Ang I]_{Tb}^{Fl}','[Ang II]_{Tb}^{Fl}','[AT1R-bound Ang II]_{Tb}^{Memb}',...
              '[AT1R-bound Ang II]_{Tb}^{Cell}','[Ang II]_{Tb}^{Cell}',...
              '[AT1R]_{Tb}^{Memb}','[AT1R]_{Tb}^{Cell}'};

PvRowNames = {'[Ang I]_{Pv}','[Ang II]_{Pv}','[AT1R-bound Ang II]_{Pv}^{Memb}'...
              '[AT1R]_{Pv}^{Memb}'};

WKRowNames = {'[Ang I]_{T}','[Ang II]_T' };

delta_SS = 100*sens_SS;
test = delta_SS;
[r,c] = find(abs(delta_SS) <= 1);
for i = 1:length(r)
    test(r(i),c(i)) = NaN;
end

test_circ = [test(1:6,:); test(8,:)]; 
test_Gl = test(9:15,:);
test_Pt = test(16:22,:);
test_Tb = test(23:29,:);
test_Pv = test(30:33,:);
test_WK = test(34:35,:);

%%
figure(10)
subplot(6,2,1)
heatmap(params(1:15),circRowNames,test_circ(:,1:15),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
    'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],...
    'ColorbarVisible','off');
ylabel('Systemic')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,2,2)
heatmap(params(16:end),circRowNames,test_circ(:,16:end),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
    'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')], 'ColorbarVisible','off');    
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

subplot(6,2,3)
heatmap(params(1:15),GlRowNames,test_Gl(:,1:15),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],'ColorbarVisible','off',...
        'ColorbarVisible','off');
ylabel('Glomerular')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,2,4)
heatmap(params(16:end),GlRowNames,test_Gl(:,16:end),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],'ColorbarVisible','off');
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

subplot(6,2,5)
heatmap(params(1:15),PtRowNames,test_Pt(:,1:15),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...%,...
       'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],...
    'ColorbarVisible','off');
ylabel('Peritubular')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,2,6)
heatmap(params(16:end),PtRowNames,test_Pt(:,16:end),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')]);
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

subplot(6,2,7)
heatmap(params(1:15),TbRowNames,test_Tb(:,1:15),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
       'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],...
    'ColorbarVisible','off');
ylabel('Tubular')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,2,8)
heatmap(params(16:end),TbRowNames,test_Tb(:,16:end),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],'ColorbarVisible','off');
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

subplot(6,2,9)
heatmap(params(1:15),PvRowNames,test_Pv(:,1:15),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...,...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],...
     'ColorbarVisible','off');  
ylabel({'Renal (blood)','vasculature'})
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,2,10)
heatmap(params(16:end),PvRowNames,test_Pv(:,16:end),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')]);
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

subplot(6,2,11)
heatmap(params(1:15),WKRowNames,test_WK(:,1:15),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...,...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],...
    'ColorbarVisible','off');  
ylabel({'Whole', 'kidney'})
xlabel('Baseline Parameters');
set(gca,'fontsize',12);

subplot(6,2,12)
heatmap(params(16:end),WKRowNames,test_WK(:,16:end),'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test(1:35,:),[],'all'),max(test(1:35,:),[],'all')],'ColorbarVisible','off');
set(gca,'fontsize',12);
xlabel('Feedback Parameters');
Ax = gca;
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

end