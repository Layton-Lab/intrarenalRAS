% Inputs:
    % days: days to run simulation
    % x0: initial condition
    % D: dose
    % units: units of dose (ng/kg/min, ng/min)
    % inf_type: type of infusion (none, SC, IV)
    % separate: whether to separate endogenous and exogenous Ang II (true, false)
function [S,t,YP] = run_model(days,x0,x_p0, D,units,inf_type,separate,varargin)
    
    % Time resolution (mins/day)
    t_res = 1440;

    % Number of points for plotting resolution
    N = ((days)*t_res)*2;
    
    % Initial time (min); Final time (min); Time vector
    t0 = 0*t_res; tend = days*t_res; tspan = linspace(t0,tend,N);
    
    % Load parameters
    pars = get_pars;
    
    % Get Dose
    if strcmp(units,'ng/kg/min')
        infusion.dose = D*955.84*0.284;
    elseif strcmp(units,'ng/min')
        infusion.dose = D*955.84;
    end
    
    infusion.type = inf_type;
    infusion.separate = separate;
    
    varargin_input = {};
    if ~isempty(varargin)
    varargin_input = varargin{:};
    end

    num_vars = 0;
    for i = 1:length(varargin_input)
        if strcmp(varargin_input{i},'losartan')
            num_vars = 26;
            break;
        end
    end
    
    if separate
        num_vars = num_vars + 75;
    else
        num_vars = num_vars + 41;
    end
    
    options = odeset('RelTol',1e-6,'AbsTol',1e-4);
    
    sol = ode15i(@(t,x,x_p) ...
                model(t,x,x_p,pars,infusion,varargin_input{:}), ...
                            tspan, x0, x_p0,options);  

    [Y,YP] = deval(sol,tspan);
    S = Y;
    t = tspan/t_res;

end