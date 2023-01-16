function p = get_pars

%%% --- Baseline (fitted) parameters

k_int         = 0.193;       % /min
k_rec         = (4/3)*0.208; % /min
k_lys         = 0.208;       % /min
k_meg         = 4.688;       % /min
k_trans       = 32.89556485; % /min
k_diff        = 0.4783;      % mL/min per g kidney
AT1R_Gl_tot   = 115;         % fmol/g kidney
c_chym        = 0.115;       % /min
c_NEP         = 0.293;       % /min
v_max         = 99.7;        % /min
k_AngI_Pt     = 6746;        % fmol/mL/min
c_ACE_Pt      = 2.0;         % /min
k_AngI_Tb     = 9992;        % fmol/mL/min
c_ACE_Tb      = 1.83;        % /min
k_AGT         = 1737;

%%% --- Parameters from the literature 

% General
W_K           = 1.49;
V_circ        = 10.3;

% Renal Volumes
V_Gl_Isf      = 0.0019;
V_Gl_cell     = 0.0019;
V_Pt_Isf      = 0.0236;
V_Tb_Pt_cell  = 0.294;
V_Tb_Fl       = 0.102;
V_Pv          = 0.085;

% Renal hemodynamics
phi_RPF       = 11.6;
FF            = 0.26;

% AT1R binding kinetics
K_D           = 1000;
k_ass         = 2.4*10^(-5);

% Systemic
h_AGT         = 240;
h_renin       = 3;
h_AngI        = 0.5;
h_AngII       = 0.267;
h_Ang17       = 0.167;
R_sec         = 1;
K_M           = 2.8*10^6;

scale = 1440/240;
k_a           = 5.41/scale/240; % /min
B_AT1R        = 2.9;
K_circ_ACE    = 3.9; % fmol/mL/min
K_circ_AGT    = 27/scale; % fmol/mL/min
K_Pt          = 4.95/scale; % fmol/mL/min
K_Tb          = 3/scale; % fmol/mL/min

AT1R_circ_tot = 5.215*10^5;
AT1R_Pv_tot   = 500;

SS = load('model_SS.mat').SSdata;
AT1R_AngII_memb_Gl_eq   = SS(11);
AT1R_AngII_memb_Pt_eq   = SS(18);
AT1R_AngII_memb_Tb_eq   = SS(25);
AT1R_AngII_memb_circ_eq = SS(4);
AngII_cell_Tb_eq        = SS(27);

p = [W_K,V_circ,V_Gl_Isf,V_Gl_cell,V_Pt_Isf,V_Tb_Pt_cell,V_Tb_Fl,V_Pv,...
     phi_RPF,FF,K_D,k_ass,h_AGT,h_renin,h_AngI,h_AngII,h_Ang17,R_sec,K_M,...
     k_int,k_rec,k_lys, k_meg, k_trans, k_diff,AT1R_Gl_tot,c_chym,c_NEP,v_max,k_AngI_Pt,...
     c_ACE_Pt,k_AngI_Tb,c_ACE_Tb,k_AGT,k_a,B_AT1R,K_circ_ACE,...
     K_circ_AGT,K_Pt,K_Tb,AT1R_circ_tot,AT1R_Pv_tot,...
     AT1R_AngII_memb_Gl_eq,AT1R_AngII_memb_Pt_eq,AT1R_AngII_memb_Tb_eq,...
     AT1R_AngII_memb_circ_eq, AngII_cell_Tb_eq];
end

