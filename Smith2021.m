
%% --- (Optional) User input --- %%

% (1) To generate select figures only, alter fig array below
    % Options:
        % 4: AGT, PRA, Ang I time series
        % 5: Ang II dose response
        % 6: Intrarenal Ang II fitting, comparing mechanism (i) and (ii)
        % 7: Intrarenal Ang II validation, comparing mechanism (i) and (ii)
        % 8: Intarenal Ang II distribution
        % 9: Steady state sensitivity analysis
        % 10: Ang II infusion sensitivity analysis
fig = [4, 5, 6, 7, 8, 9, 10]; 

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
    if fig(i) == 4 % AGT, PRA, Ang I time series
        get_Fig4_R1
    elseif fig(i) == 5 % Ang II dose response
        get_Fig5_R1
    elseif fig(i) == 6 % Intrarenal Ang II fitting, mechanism (i) vs. (ii)
        get_Fig6_R1
    elseif fig(i) == 7 % Intrarenal Ang II validation, mechanism (i) vs. (ii)
        get_Fig7_R1
    elseif fig(i) == 8 % Intrarenal Ang II distribution 
        get_Fig8_R1(dose,units,type,days)
    elseif fig(i) == 9 % Sensitivity of model steady state
        [vars_SS,sens_SS] = get_Fig9_R1(delta);
    elseif fig(i) == 10 % Senitivity of Ang II infusion results
        [vars_inf,sens_inf] = get_Fig10_R1(delta,days,dose,units,type);
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