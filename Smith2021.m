
%% --- (Optional) User input --- %%

% (1) To generate select figures only, alter fig array below
    % Options:
        % 3: AGT, PRA, Ang I time series
        % 4: Ang II dose response
        % 5: Intrarenal Ang II fitting and validation
        % 6: Intarenal Ang II distribution
        % 7: Steady state sensitivity analysis
        % 8: Ang II infusion sensitivity analysis
fig = [3,4,5,6,7,8];

% (2) To alter the dose, type, and lenngth of the Ang II experiment 
    % simulated in Figs 6 and 8, change the parameters below
    % Options:
        % dose: [0,infinity)
        % units: 'ng/min','ng/min/kg'(284 g rat)
        % type: 'SC','IV'
        % days: [0,infinity)
dose = 40; units = 'ng/min'; type = 'SC'; days = 13; % default

% (3) To alter the percentage used in the sensitivity analysis
    % change delta ([-1,1) accordingly
delta = 0.1; % 10 percent increase (default)
        
%% --- Generate Figures --- %%

for i = 1:length(fig)
    if fig(i) == 3 % AGT, PRA, Ang I time series
        get_Fig3
    elseif fig(i) == 4 % Ang II dose response
        get_Fig4
    elseif fig(i) == 5 % Intrarenal Ang II fitting and validation
        get_Fig5
    elseif fig(i) == 6 % Intrarenal Ang II distribution 
        get_Fig6(dose,units,type,days)
    elseif fig(i) == 7 % Sensitivity of model steady state
        [vars_SS,sens_SS] = get_Fig7(delta);
    elseif fig(i) == 8 % Senitivity of Ang II infusion results
        [vars_inf,sens_inf] = get_Fig8(delta,days,dose,units,type);
    end
end

%% --- Save figures as PDFs --- %%

% (4) To save desired figure as pdf
    % (i) Click on figure you want to save
    % (ii) Set save_fig = true
    % (iii) run the folowing with desired title as input  
title = 'Fig';
save_fig = false;
if save_fig
    save_pdf(title)
end