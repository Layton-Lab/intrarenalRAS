%% --- (Optional) User input --- %%

% (1) To generate select figures only, alter fig array below
    % Options:
        % "2": Model fit to a single oral dose of Losartan
        % "3-5": Replicating Losartan experiments from Zou et al.
          % 3: Systemic and intrarenal RAS peptide changes 
          % 4: Plasma and intrarenal endo vs. exo Ang II accumulation
          % 5: Intrarenal distribution of endo vs. exo Ang II
        % "6": Effect of Losartan treatment post hypetension-induction
        % "7-8": sensitivity analysis
          % 7: Effect of parameter perturbations on Losartan and EXP3174
          % 8: Effect of parameter perturbations on the RAS
fig = ["2", "3-5", "6", "7-8"];

% (2) To alter the dose, type, and length of the Ang II/Losartan experiment 
    % simulated in Figs 6-8, change the parameters below
    % Options:
        % dose: [0,infinity)
        % losartan_dose: [0,infinity)
        % units: 'ng/min','ng/min/kg'(284 g rat)
        % type: 'SC','IV'
        % days: [0,infinity)
        % days: [1, infinity)
dose = 40; units = 'ng/min'; type = 'SC'; days = 13; % default
losartan_dose = 30; days_losartan = 7; % default

% (3) To alter the percentage used in the sensitivity analysis
    % change delta ([-1,1]) accordingly
delta = 0.1; % 10 percent increase (default)
        
%% --- Generate Figures --- %%

for i = 1:length(fig)
    if strcmp(fig(i),"2")
        get_Losartan_Fig2
    elseif strcmp(fig(i),"3-5") 
        get_Losartan_Fig3_to_5
    elseif strcmp(fig(i), "6") 
        get_Losartan_Fig6(days,dose,units,type,losartan_dose, days_losartan)
    elseif strcmp(fig(i),"7-8") 
        get_Losartan_Fig7_and_8(delta,days,dose,units,type,losartan_dose)
    else
        disp('Unknown input.')
    end
end

%% --- Save figures as PDFs --- %%

% (4) To save desired figure as pdf
    % (i) Click on figure you want to save
    % (ii) Set save_fig = true
    % (iii) run the folowing with desired title as input  
tit = 'Figure';
save_fig = false;
if save_fig
    save_pdf(tit)
end