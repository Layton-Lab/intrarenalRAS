
% Fig 7: sensitivity of the model steady state
function [SSdata,test] = get_Fig9_R1(delta)

% Don't infuse Ang II
infusion.dose = 0; infusion.type = 'none';

% Don't separate exogenous from endogenous Ang II
infusion.separate = false;

% Turn off all feedback
varargin_input = {'all','renin'}; 

% Number of days to run the simulation; time resolution (number of minutes/day)
days = 2; t_res = 1440; 
    
% Number of points for plotting resolution
N = ((days)*t_res);
    
% Initial time (min); Final time (min); Time vector
t0 = 0*t_res; tend = days*t_res; tspan = linspace(t0,tend,N);

% Load original parameters
pars = get_pars;
    
% Load initial condition
x0 = [load('model_SS.mat').SSdata];
num_vars = length(x0); x_p0 = zeros(num_vars,1);

% Set solver options
options = odeset('RelTol',1e-6,'AbsTol',1e-4);

% Number of baseline parameters 
p_index = 15;

% Initialize ,matrix to hold new SS for each parameter change
SSdata = zeros(num_vars,p_index);

% Initialize matrix to hold % change in each SS for each parameter change
sens_SS = zeros(num_vars,p_index);

for i = 1:p_index  
     % Increase baseline parameter i by delta %
     pars_i = pars;
     delta_pi = delta*pars(i+19); % Baseline parameters start at index 20 in pars
     pars_i(i+19) = pars(i+19) + delta_pi;
     
     % Solve new SS
     sol_0 = ode15i(@(t,x,x_p) ...
                model(t,x,x_p,pars_i,infusion,varargin_input{:}), ...
                            tspan, x0, x_p0,options);  
     sim = sol_0.y;
     
     % Store new SS and % change in SS
     SSdata(:,i) = sim(:,end);
     sens_SS(:,i) = (SSdata(:,i) - x0)./x0;
end

%% Plot table

% Parameter labels
params = {'k_{int}','k_{rec}', 'k_{lys}', 'k_{meg}', 'k_{sec}', 'k_{diff}', 'AT1R_{Gl}^{tot}','c_{chym}','c_{NEP}',...
          'v_{max}','k_{AngI_{Pt}}','c_{ACE_{Pt}}','k_{AngI_{Tb}}','c_{ACE_{Tb}}','k_{AGT}'};

% Variable labels
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

% Set all values <1% in magnitude to NaN for plotting
delta_SS = 100*sens_SS;
test = delta_SS;
[r,c] = find(abs(delta_SS) <= 1);
for i = 1:length(r)
    test(r(i),c(i)) = NaN;
end

% Separate each compartment
test_circ = [test(1:6,:); test(8,:)];
test_Gl = test(9:15,:);
test_Pt = test(16:22,:);
test_Tb = test(23:29,:);
test_Pv = test(30:33,:);
test_WK = test(34:35,:);

figure(9)
subplot(6,1,1)
heatmap(params,circRowNames,test_circ,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
    'ColorLimits',[min(test,[],'all'),max(test,[],'all')],'colorbarvisible','off');
ylabel('Systemic')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,1,2)
heatmap(params,GlRowNames,test_Gl,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
        'ColorLimits',[min(test,[],'all'),max(test,[],'all')],'colorbarvisible','off');
ylabel('Glomerular')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));


subplot(6,1,3)
heatmap(params,PtRowNames,test_Pt,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
       'ColorLimits',[min(test,[],'all'),max(test,[],'all')]);
ylabel('Peritubular')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,1,4)
heatmap(params,TbRowNames,test_Tb,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...
       'ColorLimits',[min(test,[],'all'),max(test,[],'all')],'colorbarvisible','off');
ylabel('Tubular')
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

subplot(6,1,5)
heatmap(params,PvRowNames,test_Pv,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...,...
        'ColorLimits',[min(test,[],'all'),max(test,[],'all')],'colorbarvisible','off');
ylabel({'Renal (blood)','vasculature'})
set(gca,'fontsize',12);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));


subplot(6,1,6)
heatmap(params,WKRowNames,test_WK,'Colormap',summer,...
    'MissingDataColor','w','MissingDataLabel','<1%',...,...
        'ColorLimits',[min(test,[],'all'),max(test,[],'all')],'colorbarvisible','off');
ylabel({'Whole', 'kidney'})
xlabel('Parameters');
set(gca,'fontsize',12);

end