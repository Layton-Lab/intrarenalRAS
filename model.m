function f = model(t,x,x_p,pars,infusion,varargin) 
%% Parameters

lit_pars      = pars(1:19);
base_pars     = pars(19+1:19+1+13);
fb_pars       = pars(19+1+13:19+1+13+8);
SS            = pars(19+1+13+8:end);

%%% --- Baseline (fitted) parameters
k_int         = base_pars(1);
k_rec         = base_pars(2); 
k_lys         = base_pars(3);
AT1R_Gl_tot   = base_pars(4);
c_chym        = base_pars(5); c_ACE_circ = 9*c_chym;
c_NEP         = base_pars(6); c_ACE2     = 0.21*c_NEP;
v_max         = base_pars(7);
k_AngI_Pt     = base_pars(8);
c_ACE_Pt      = base_pars(9);
k_AngI_Tb     = base_pars(10);
c_ACE_Tb      = base_pars(11);
S_Tb          = base_pars(12);
k_AGT         = base_pars(13);

%%% --- Parameters from the literature 

% General
W_K           = lit_pars(1);
V_circ        = lit_pars(2);

% Renal Volumes
V_Gl_Isf      = lit_pars(3);
V_Gl_cell     = lit_pars(4);
V_Pt_Isf      = lit_pars(5);
V_Tb_Pt_cell  = lit_pars(6);
V_Tb_Fl       = lit_pars(7);
V_Pv          = lit_pars(8);

% Renal hemodynamics
phi_RPF       = lit_pars(9)/W_K;
FF            = lit_pars(10);
phi_GFR       = FF*phi_RPF;
phi_L         = (0.02/1.02)*phi_GFR; 
phi_U         = phi_L; 
phi_Pt        = S_Tb*(phi_GFR-phi_U);
phi_Pv        = phi_GFR - phi_U; 

% AT1R binding kinetics
K_D           = lit_pars(11);
k_ass         = lit_pars(12);
k_diss        = k_ass*K_D;

% Systemic
h_AGT         = lit_pars(13); v_AGT   = log(2)/h_AGT;
h_renin       = lit_pars(14); v_renin = log(2)/h_renin;
h_AngI        = lit_pars(15); v_I     = log(2)/h_AngI;
h_AngII       = lit_pars(16); v_II    = log(2)/h_AngII;
h_Ang17       = lit_pars(17); v_17    = log(2)/h_Ang17;
R_sec         = lit_pars(18);
K_M           = lit_pars(19);

% Feedback parameters
k_a           = fb_pars(1)*10^(-3);
B_AT1R        = fb_pars(2);
K_circ_ACE    = fb_pars(3);
K_circ_AGT    = fb_pars(4)*10^2;
K_Pt          = fb_pars(5);
K_Tb          = fb_pars(6)*10^(-1);
AT1R_circ_tot = fb_pars(7);
AT1R_Pv_tot   = fb_pars(8);

% Steady states 
AT1R_AngII_memb_Gl_eq   = SS(1);
AT1R_AngII_memb_Pt_eq   = SS(2);
AT1R_AngII_memb_Tb_eq   = SS(3);
AT1R_AngII_memb_circ_eq = SS(4);

%% Optional inputs
setNuTo1 = false;
for i = 1:length(varargin)
    if strcmp(varargin{i},'all') % remove all fb
        K_circ_ACE    = 0;
        K_circ_AGT    = 0;
        K_Pt          = 0;
        K_Tb          = 0;
    elseif strcmp(varargin{i},'all_renal') % remove all intrarenal fb
        K_Pt = 0;
        K_Tb = 0;
    elseif strcmp(varargin{i},'circ_ACE') % remove systemic ACE fb
        K_circ_ACE    = 0;
    elseif strcmp(varargin{i},'circ_AGT') % remove hepatic AGT fb
        K_circ_AGT    = 0;
    elseif strcmp(varargin{i},'Pt') % remove basolateral AT1R expression fb
        K_Pt = 0;
    elseif strcmp(varargin{i},'Tb') % remove apical AT1R expression fb
        K_Tb = 0;
    elseif strcmp(varargin{i},'renin') % remove renin secretion fb
        setNuTo1 = true;
    elseif strcmp(varargin{i},'new_SS') % for infusion sensitivity analysis
        AT1R_AngII_memb_Gl_eq = varargin{i+1};
        AT1R_AngII_memb_Pt_eq   = varargin{i+2};
        AT1R_AngII_memb_Tb_eq   = varargin{i+3};
        AT1R_AngII_memb_circ_eq = varargin{i+4};
        i = i+1;
    end
end
%% Variables

% Systemic 
AGT_circ             = x(1); AGT_circ_p             = x_p(1);
AngI_circ            = x(2); AngI_circ_p            = x_p(2);
AngII_circ           = x(3); AngII_circ_p           = x_p(3);
AT1R_AngII_memb_circ = x(4); AT1R_AngII_memb_circ_p = x_p(4);
AT1R_memb_circ       = x(5); AT1R_memb_circ_p       = x_p(5);
Ang17_circ           = x(6); Ang17_circ_p           = x_p(6);
PRC                  = x(7); PRC_p                  = x_p(7);
PRA                  = x(8); PRA_p                  = x_p(8);

% Glomerular
AngI_Isf_Gl          = x(9); AngI_Isf_Gl_p          = x_p(9);
AngII_Isf_Gl         = x(10); AngII_Isf_Gl_p        = x_p(10);
AT1R_AngII_memb_Gl   = x(11); AT1R_AngII_memb_Gl_p  = x_p(11);
AT1R_AngII_cell_Gl   = x(12); AT1R_AngII_cell_Gl_p  = x_p(12);
AngII_cell_Gl        = x(13); AngII_cell_Gl_p       = x_p(13);
AT1R_memb_Gl         = x(14); AT1R_memb_Gl_p        = x_p(14);
AT1R_cell_Gl         = x(15); AT1R_cell_Gl_p        = x_p(15);

% Peritubular
AngI_Isf_Pt          = x(16); AngI_Isf_Pt_p         = x_p(16);
AngII_Isf_Pt         = x(17); AngII_Isf_Pt_p        = x_p(17);
AT1R_AngII_memb_Pt   = x(18); AT1R_AngII_memb_Pt_p  = x_p(18);
AT1R_AngII_cell_Pt   = x(19); AT1R_AngII_cell_Pt_p  = x_p(19);
AngII_cell_Pt        = x(20); AngII_cell_Pt_p       = x_p(20);
AT1R_memb_Pt         = x(21); AT1R_memb_Pt_p        = x_p(21);
AT1R_cell_Pt         = x(22); AT1R_cell_Pt_p        = x_p(22);

% Tubular
AngI_Fl_Tb           = x(23); AngI_Fl_Tb_p          = x_p(23);
AngII_Fl_Tb          = x(24); AngII_Fl_Tb_p         = x_p(24);
AT1R_AngII_memb_Tb   = x(25); AT1R_AngII_memb_Tb_p  = x_p(25);
AT1R_AngII_cell_Tb   = x(26); AT1R_AngII_cell_Tb_p  = x_p(26);
AngII_cell_Tb        = x(27); AngII_cell_Tb_p       = x_p(27);
AT1R_memb_Tb         = x(28); AT1R_memb_Tb_p        = x_p(28);
AT1R_cell_Tb         = x(29); AT1R_cell_Tb_p        = x_p(29);

% Renal (blood) vasculature
AngI_Pv              = x(30); AngI_Pv_p             = x_p(30);
AngII_Pv             = x(31); AngII_Pv_p            = x_p(31);
AT1R_AngII_memb_Pv   = x(32); AT1R_AngII_memb_Pv_p  = x_p(32);
AT1R_memb_Pv         = x(33); AT1R_memb_Pv_p        = x_p(33);

% Whole Kidney
AngI_T               = x(34); AngI_T_p              = x_p(34);
AngII_T              = x(35); AngII_T_p             = x_p(35);

% Feedback
nu_AT1R              = x(36); nu_AT1R_p             = x_p(36);
fb_circ_AGT          = x(37); fb_circ_AGT_p         = x_p(37);
fb_circ_ACE          = x(38); fb_circ_ACE_p         = x_p(38);
fb_Pt                = x(39); fb_Pt_p               = x_p(39);
fb_Tb                = x(40); fb_Tb_p               = x_p(40);

if infusion.separate
    
    %%% --- Exogenous 
    
    % Systemic
    AngII_circ_exo           = x(41); AngII_circ_exo_p           = x_p(41);
    AT1R_AngII_memb_circ_exo = x(42); AT1R_AngII_memb_circ_exo_p = x_p(42);
    
    % Glomerular
    AngII_Isf_Gl_exo         = x(43); AngII_Isf_Gl_exo_p         = x_p(43);
    AT1R_AngII_memb_Gl_exo   = x(44); AT1R_AngII_memb_Gl_exo_p   = x_p(44);
    AT1R_AngII_cell_Gl_exo   = x(45); AT1R_AngII_cell_Gl_exo_p   = x_p(45);
    AngII_cell_Gl_exo        = x(46); AngII_cell_Gl_exo_p        = x_p(46);
    
    % Peritubular
    AngII_Isf_Pt_exo         = x(47); AngII_Isf_Pt_exo_p         = x_p(47);
    AT1R_AngII_memb_Pt_exo   = x(48); AT1R_AngII_memb_Pt_exo_p   = x_p(48);
    AT1R_AngII_cell_Pt_exo   = x(49); AT1R_AngII_cell_Pt_exo_p   = x_p(49);
    AngII_cell_Pt_exo        = x(50); AngII_cell_Pt_exo_p        = x_p(50);
    
    % Tubular
    AngII_Fl_Tb_exo          = x(51); AngII_Fl_Tb_exo_p         = x_p(51);
    AT1R_AngII_memb_Tb_exo   = x(52); AT1R_AngII_memb_Tb_exo_p  = x_p(52);
    AT1R_AngII_cell_Tb_exo   = x(53); AT1R_AngII_cell_Tb_exo_p  = x_p(53);
    AngII_cell_Tb_exo        = x(54); AngII_cell_Tb_exo_p       = x_p(54);
    
    % Renal (blood) vasculatur
    AngII_Pv_exo             = x(55); AngII_Pv_exo_p            = x_p(55);
    AT1R_AngII_memb_Pv_exo   = x(56); AT1R_AngII_memb_Pv_exo_p  = x_p(56);

    % Whole Kidney
    AngII_T_exo              = x(57); AngII_T_exo_p             = x_p(57);
    
    %%%--- Total
    
    % Systemic
    AngII_circ_tot           = x(58); AngII_circ_tot_p           = x_p(58);
    AT1R_AngII_memb_circ_tot = x(59); AT1R_AngII_memb_circ_tot_p = x_p(59);
    
    % Glomerular
    AngII_Isf_Gl_tot         = x(60); AngII_Isf_Gl_tot_p         = x_p(60);
    AT1R_AngII_memb_Gl_tot   = x(61); AT1R_AngII_memb_Gl_tot_p   = x_p(61);
    AT1R_AngII_cell_Gl_tot   = x(62); AT1R_AngII_cell_Gl_tot_p   = x_p(62);
    AngII_cell_Gl_tot        = x(63); AngII_cell_Gl_tot_p        = x_p(63);
    
    % Peritubular
    AngII_Isf_Pt_tot         = x(64); AngII_Isf_Pt_tot_p         = x_p(64);
    AT1R_AngII_memb_Pt_tot   = x(65); AT1R_AngII_memb_Pt_tot_p   = x_p(65);
    AT1R_AngII_cell_Pt_tot   = x(66); AT1R_AngII_cell_Pt_tot_p   = x_p(66);
    AngII_cell_Pt_tot        = x(67); AngII_cell_Pt_tot_p        = x_p(67);
    
    % Tubular
    AngII_Fl_Tb_tot          = x(68); AngII_Fl_Tb_tot_p         = x_p(68);
    AT1R_AngII_memb_Tb_tot   = x(69); AT1R_AngII_memb_Tb_tot_p  = x_p(69);
    AT1R_AngII_cell_Tb_tot   = x(70); AT1R_AngII_cell_Tb_tot_p  = x_p(70);
    AngII_cell_Tb_tot        = x(71); AngII_cell_Tb_tot_p       = x_p(71);
    
    % Renal (blood) vasculature
    AngII_Pv_tot             = x(72); AngII_Pv_tot_p            = x_p(72);
    AT1R_AngII_memb_Pv_tot   = x(73); AT1R_AngII_memb_Pv_tot_p  = x_p(73);

    % Whole Kidney
    AngII_T_tot              = x(74); AngII_T_tot_p             = x_p(74);
end
%% Infusion

D = infusion.dose;
K_AngII = 0;
if strcmp(infusion.type,'SC')
    K_AngII = (D/V_circ)*(1-exp(-k_a*t));
elseif  strcmp(infusion.type,'IV')
    K_AngII = D/V_circ;
end

%% Feedback

if infusion.separate
    q_circ = AT1R_AngII_memb_circ_tot/AT1R_AngII_memb_circ_eq;
    q_Gl = AT1R_AngII_memb_Gl_tot/AT1R_AngII_memb_Gl_eq;
    q_Pt = AT1R_AngII_memb_Pt_tot/AT1R_AngII_memb_Pt_eq;
    q_Tb = AT1R_AngII_memb_Tb_tot/AT1R_AngII_memb_Tb_eq;
else
    q_circ = AT1R_AngII_memb_circ/AT1R_AngII_memb_circ_eq;
    q_Gl = AT1R_AngII_memb_Gl/AT1R_AngII_memb_Gl_eq;
    q_Pt = AT1R_AngII_memb_Pt/AT1R_AngII_memb_Pt_eq;
    q_Tb = AT1R_AngII_memb_Tb/AT1R_AngII_memb_Tb_eq;
end
%% Model

f = zeros(length(x),1);

%%% --- Systemic 

f(1) = AGT_circ_p  - (k_AGT + fb_circ_AGT - PRA - v_AGT*AGT_circ);

f(2) = AngI_circ_p - (PRA - (c_chym + c_ACE_circ + fb_circ_ACE + c_NEP + v_I)*AngI_circ...
                   + (W_K/V_circ)*phi_L*(AngI_Isf_Pt + AngI_Isf_Gl)...
                   + (W_K/V_circ)*((phi_RPF-2*phi_L-phi_U)*AngI_Pv - phi_RPF*AngI_circ));

if infusion.separate
    f(3) = AngII_circ_p - ((c_chym + c_ACE_circ + fb_circ_ACE)*AngI_circ...  
                        - (c_ACE2 + v_II)*AngII_circ...
                        + (W_K/V_circ)*phi_L*(AngII_Isf_Pt+AngII_Isf_Gl)...
                        + (W_K/V_circ)*(phi_RPF-2*phi_L-phi_U)*AngII_Pv ...
                        - (W_K/V_circ)*phi_RPF*AngII_circ...
                        + k_diss*AT1R_AngII_memb_circ - k_ass*AT1R_memb_circ*AngII_circ);              
else
    f(3) = AngII_circ_p - (K_AngII + (c_chym + c_ACE_circ + fb_circ_ACE)*AngI_circ...  
                        - (c_ACE2 + v_II)*AngII_circ...
                        + (W_K/V_circ)*phi_L*(AngII_Isf_Pt+AngII_Isf_Gl)...
                        + (W_K/V_circ)*(phi_RPF-2*phi_L-phi_U)*AngII_Pv ...
                        - (W_K/V_circ)*phi_RPF*AngII_circ...
                        + k_diss*AT1R_AngII_memb_circ - k_ass*AT1R_memb_circ*AngII_circ); 
end

f(4) = AT1R_AngII_memb_circ_p - (k_ass*AT1R_memb_circ*AngII_circ ...
                              - k_diss*AT1R_AngII_memb_circ);
                          
if infusion.separate
    f(5) = AT1R_memb_circ - (AT1R_circ_tot - AT1R_AngII_memb_circ_tot);
else
    f(5) = AT1R_memb_circ - (AT1R_circ_tot - AT1R_AngII_memb_circ);
end

f(6) = Ang17_circ_p - (c_NEP*AngI_circ + c_ACE2*AngII_circ - v_17*Ang17_circ);

f(7) = PRC_p - (R_sec*nu_AT1R - v_renin*PRC);

f(8) = PRA - v_max*PRC*AGT_circ/(AGT_circ+K_M);

%%% ----  Glomerular

f(9)  = AngI_Isf_Gl_p  - ((phi_L/V_Gl_Isf)*(AngI_circ -AngI_Isf_Gl));

f(10) = AngII_Isf_Gl_p - ((phi_L/V_Gl_Isf)*(AngII_circ - AngII_Isf_Gl) ...
                       + k_diss*AT1R_AngII_memb_Gl -k_ass*AngII_Isf_Gl*AT1R_memb_Gl);
                   
f(11) = AT1R_AngII_memb_Gl_p - (k_ass*AngII_Isf_Gl*AT1R_memb_Gl ...
                             - (k_diss + k_int)*AT1R_AngII_memb_Gl);

f(12) = AT1R_AngII_cell_Gl_p - ((V_Gl_Isf/V_Gl_cell)*k_int*AT1R_AngII_memb_Gl ...
                             + k_ass*AngII_cell_Gl*AT1R_cell_Gl - k_diss*AT1R_AngII_cell_Gl) ;

f(13) = AngII_cell_Gl_p - (k_diss*AT1R_AngII_cell_Gl - k_ass*AngII_cell_Gl*AT1R_cell_Gl ...
                        - k_lys*AngII_cell_Gl);
                    
if infusion.separate                   
    f(14) = AT1R_memb_Gl - (AT1R_Gl_tot/V_Gl_Isf - AT1R_AngII_memb_Gl_tot ...
                         - (V_Gl_cell/V_Gl_Isf)*(AT1R_AngII_cell_Gl_tot + AT1R_cell_Gl));
else
    f(14) = AT1R_memb_Gl - (AT1R_Gl_tot/V_Gl_Isf - AT1R_AngII_memb_Gl ...
                         - (V_Gl_cell/V_Gl_Isf)*(AT1R_AngII_cell_Gl + AT1R_cell_Gl));    
end

if infusion.separate
    f(15) = AT1R_cell_Gl_p - (k_diss*AT1R_AngII_cell_Gl_tot - k_ass*AngII_cell_Gl_tot*AT1R_cell_Gl ...
                           - k_rec*AT1R_cell_Gl);
else
    f(15) = AT1R_cell_Gl_p - (k_diss*AT1R_AngII_cell_Gl - k_ass*AngII_cell_Gl*AT1R_cell_Gl ...
                           - k_rec*AT1R_cell_Gl);   
end

%%% ----  Peritubular

f(16) = AngI_Isf_Pt_p  - (k_AngI_Pt + (phi_Pt/V_Pt_Isf)*AngI_Fl_Tb...
                       - (c_ACE_Pt + (phi_Pv+phi_L)/V_Pt_Isf)*AngI_Isf_Pt);
                    
f(17) = AngII_Isf_Pt_p - (c_ACE_Pt*AngI_Isf_Pt + (phi_Pt/V_Pt_Isf)*AngII_Fl_Tb ...
                       - ((phi_Pv+phi_L)/V_Pt_Isf)*AngII_Isf_Pt ...
                       + k_diss*AT1R_AngII_memb_Pt - k_ass*AngII_Isf_Pt*AT1R_memb_Pt);
                   
f(18) = AT1R_AngII_memb_Pt_p - (k_ass*AngII_Isf_Pt*AT1R_memb_Pt ...
                             - (k_diss + k_int)*AT1R_AngII_memb_Pt);

f(19) = AT1R_AngII_cell_Pt_p - ((V_Pt_Isf/V_Tb_Pt_cell)*k_int*AT1R_AngII_memb_Pt ...
                             + k_ass*AngII_cell_Pt*AT1R_cell_Pt - k_diss*AT1R_AngII_cell_Pt) ;

f(20) = AngII_cell_Pt_p - (k_diss*AT1R_AngII_cell_Pt - k_ass*AngII_cell_Pt*AT1R_cell_Pt ...
                        - k_lys*AngII_cell_Pt);
                    
if infusion.separate                   
    f(21) = AT1R_memb_Pt_p - (fb_Pt + (V_Tb_Pt_cell/V_Pt_Isf)*k_rec*AT1R_cell_Pt ...
                           + k_diss*AT1R_AngII_memb_Pt_tot ...
                           - k_ass*AngII_Isf_Pt_tot*AT1R_memb_Pt);           
else
    f(21) = AT1R_memb_Pt_p - (fb_Pt + (V_Tb_Pt_cell/V_Pt_Isf)*k_rec*AT1R_cell_Pt ...
                           + k_diss*AT1R_AngII_memb_Pt ...
                           - k_ass*AngII_Isf_Pt*AT1R_memb_Pt);   
end

if infusion.separate
    f(22) = AT1R_cell_Pt_p - (k_diss*AT1R_AngII_cell_Pt_tot - k_ass*AngII_cell_Pt_tot*AT1R_cell_Pt ...
                           - k_rec*AT1R_cell_Pt);
else
    f(22) = AT1R_cell_Pt_p - (k_diss*AT1R_AngII_cell_Pt - k_ass*AngII_cell_Pt*AT1R_cell_Pt ...
                           - k_rec*AT1R_cell_Pt);   
end

%%% ----  Tubular

f(23) = AngI_Fl_Tb_p  - (k_AngI_Tb + (phi_GFR/V_Tb_Fl)*AngI_circ ...
                      - (c_ACE_Tb + (phi_U + phi_Pt)/V_Tb_Fl)*AngI_Fl_Tb);
                    
f(24) = AngII_Fl_Tb_p - ((phi_GFR/V_Tb_Fl)*AngII_circ + c_ACE_Tb*AngI_Fl_Tb ...
                      + k_diss*AT1R_AngII_memb_Tb - k_ass*AngII_Fl_Tb*AT1R_memb_Tb ...
                      - ((phi_U + phi_Pt)/V_Tb_Fl)*AngII_Fl_Tb);

f(25) = AT1R_AngII_memb_Tb_p - (k_ass*AngII_Fl_Tb*AT1R_memb_Tb ...
                             - (k_diss + k_int)*AT1R_AngII_memb_Tb);

f(26) = AT1R_AngII_cell_Tb_p - ((V_Tb_Fl/V_Tb_Pt_cell)*k_int*AT1R_AngII_memb_Tb ...
                             + k_ass*AngII_cell_Tb*AT1R_cell_Tb - k_diss*AT1R_AngII_cell_Tb) ;

f(27) = AngII_cell_Tb_p - (k_diss*AT1R_AngII_cell_Tb - k_ass*AngII_cell_Tb*AT1R_cell_Tb ...
                        - k_lys*AngII_cell_Tb);
                    
if infusion.separate                  
    f(28) = AT1R_memb_Tb_p - (fb_Tb + (V_Tb_Pt_cell/V_Tb_Fl)*k_rec*AT1R_cell_Tb ...
                           + k_diss*AT1R_AngII_memb_Tb_tot ...
                           - k_ass*AngII_Fl_Tb_tot*AT1R_memb_Tb); 
else
    f(28) = AT1R_memb_Tb_p - (fb_Tb + (V_Tb_Pt_cell/V_Tb_Fl)*k_rec*AT1R_cell_Tb ...
                           + k_diss*AT1R_AngII_memb_Tb ...
                           - k_ass*AngII_Fl_Tb*AT1R_memb_Tb);   
end

if infusion.separate
    f(29) = AT1R_cell_Tb_p - (k_diss*AT1R_AngII_cell_Tb_tot - k_ass*AngII_cell_Tb_tot*AT1R_cell_Tb ...
                           - k_rec*AT1R_cell_Tb);
else
    f(29) = AT1R_cell_Tb_p - (k_diss*AT1R_AngII_cell_Tb - k_ass*AngII_cell_Tb*AT1R_cell_Tb ...
                           - k_rec*AT1R_cell_Tb);   
end

%%% --- Renal (blood) vasculature

f(30) = AngI_Pv_p - (((phi_RPF-phi_GFR-phi_L)/V_Pv)*AngI_circ ...
                  + (phi_Pv/V_Pv)*AngI_Isf_Pt ... 
                  - (((phi_RPF-2*phi_L-phi_U)/V_Pv) + v_I)*AngI_Pv);

f(31) = AngII_Pv_p - (((phi_RPF - phi_GFR - phi_L)/V_Pv)*AngII_circ ...
                   + (phi_Pv/V_Pv)*AngII_Isf_Pt...
                   - (((phi_RPF-2*phi_L-phi_U)/V_Pv)+v_II)*AngII_Pv...
                   - k_ass*AT1R_memb_Pv*AngII_Pv + k_diss*AT1R_AngII_memb_Pv); 
               
f(32) = AT1R_AngII_memb_Pv_p - (k_ass*AT1R_memb_Pv*AngII_Pv - k_diss*AT1R_AngII_memb_Pv);

if infusion.separate 
    f(33) = AT1R_memb_Pv - (AT1R_Pv_tot - AT1R_AngII_memb_Pv_tot);
else
    f(33) = AT1R_memb_Pv - (AT1R_Pv_tot - AT1R_AngII_memb_Pv);
end

%%% --- Whole Kidney

f(34) = AngI_T - (V_Gl_Isf*AngI_Isf_Gl + V_Pt_Isf*AngI_Isf_Pt ...
               + V_Tb_Fl*AngI_Fl_Tb + V_Pv*AngI_Pv);
           
f(35) = AngII_T - (V_Gl_Isf*AngII_Isf_Gl + V_Pt_Isf*AngII_Isf_Pt + V_Tb_Fl*AngII_Fl_Tb...
                + V_Gl_Isf*AT1R_AngII_memb_Gl + V_Pt_Isf*AT1R_AngII_memb_Pt + V_Tb_Fl*AT1R_AngII_memb_Tb...
                + V_Tb_Pt_cell*(AT1R_AngII_cell_Tb + AT1R_AngII_cell_Pt + AngII_cell_Tb + AngII_cell_Pt)...
                + V_Gl_cell*(AT1R_AngII_cell_Gl + AngII_cell_Gl) + V_Pv*(AngII_Pv+AT1R_AngII_memb_Pv));
                     
%%% --- Feedback

if setNuTo1 
    f(36) = nu_AT1R - 1;
else
    f(36) = nu_AT1R - q_Gl^(-B_AT1R);
end

if q_circ <= 1 
    f(37) = fb_circ_AGT - 0;
    f(38) = fb_circ_ACE - 0;
elseif q_circ > 1
    f(37) = fb_circ_AGT - (K_circ_AGT*(q_circ-1));
    f(38) = fb_circ_ACE - (K_circ_ACE*(q_circ-1));
end

if  q_Pt <= 1 
    f(39) = fb_Pt - 0;
elseif q_Pt > 1
    f(39) = fb_Pt - (K_Pt*(q_Pt-1));
end

if q_Tb <= 1  
    f(40) = fb_Tb - 0;
elseif q_Pt > 1
    f(40) = fb_Tb - (K_Tb*(q_Tb-1));
end


if infusion.separate
    
    %%% --- Exogenous
    
    % Systemic
    f(41) = AngII_circ_exo_p - (K_AngII - (c_ACE2 + v_II)*AngII_circ_exo ...
                             + (W_K/V_circ)*phi_L*(AngII_Isf_Pt_exo+AngII_Isf_Gl_exo)...
                             + (W_K/V_circ)*(phi_RPF-2*phi_L-phi_U)*AngII_Pv_exo ...
                             - (W_K/V_circ)*phi_RPF*AngII_circ_exo...
                             + k_diss*AT1R_AngII_memb_circ_exo - k_ass*AT1R_memb_circ*AngII_circ_exo); 
                         
    f(42) = AT1R_AngII_memb_circ_exo_p - (k_ass*AT1R_memb_circ*AngII_circ_exo ...
                                       - k_diss*AT1R_AngII_memb_circ_exo);
    
    % Glomerular
    f(43) = AngII_Isf_Gl_exo_p - ((phi_L/V_Gl_Isf)*(AngII_circ_exo - AngII_Isf_Gl_exo) ...
                       + k_diss*AT1R_AngII_memb_Gl_exo -k_ass*AngII_Isf_Gl_exo*AT1R_memb_Gl);   
                   
    f(44) = AT1R_AngII_memb_Gl_exo_p - (k_ass*AngII_Isf_Gl_exo*AT1R_memb_Gl ...
                                     - (k_diss + k_int)*AT1R_AngII_memb_Gl_exo);

    f(45) = AT1R_AngII_cell_Gl_exo_p - ((V_Gl_Isf/V_Gl_cell)*k_int*AT1R_AngII_memb_Gl_exo ...
                                     + k_ass*AngII_cell_Gl_exo*AT1R_cell_Gl - k_diss*AT1R_AngII_cell_Gl_exo) ;

    f(46) = AngII_cell_Gl_exo_p - (k_diss*AT1R_AngII_cell_Gl_exo - k_ass*AngII_cell_Gl_exo*AT1R_cell_Gl ...
                                - k_lys*AngII_cell_Gl_exo);
    
    % Peritubular 
    f(47) = AngII_Isf_Pt_exo_p - ((phi_Pt/V_Pt_Isf)*AngII_Fl_Tb_exo ...
                               - ((phi_Pv+phi_L)/V_Pt_Isf)*AngII_Isf_Pt_exo ...
                               + k_diss*AT1R_AngII_memb_Pt_exo - k_ass*AngII_Isf_Pt_exo*AT1R_memb_Pt);
                   
    f(48) = AT1R_AngII_memb_Pt_exo_p - (k_ass*AngII_Isf_Pt_exo*AT1R_memb_Pt ...
                                     - (k_diss + k_int)*AT1R_AngII_memb_Pt_exo);

    f(49) = AT1R_AngII_cell_Pt_exo_p - ((V_Pt_Isf/V_Tb_Pt_cell)*k_int*AT1R_AngII_memb_Pt_exo ...
                                     + k_ass*AngII_cell_Pt_exo*AT1R_cell_Pt - k_diss*AT1R_AngII_cell_Pt_exo) ;

    f(50) = AngII_cell_Pt_exo_p - (k_diss*AT1R_AngII_cell_Pt_exo - k_ass*AngII_cell_Pt_exo*AT1R_cell_Pt ...
                                - k_lys*AngII_cell_Pt_exo);
    
    % Tubular
    f(51) = AngII_Fl_Tb_exo_p - ((phi_GFR/V_Tb_Fl)*AngII_circ_exo ...
                              + k_diss*AT1R_AngII_memb_Tb_exo - k_ass*AngII_Fl_Tb_exo*AT1R_memb_Tb ...
                              - ((phi_U + phi_Pt)/V_Tb_Fl)*AngII_Fl_Tb_exo);

    f(52) = AT1R_AngII_memb_Tb_exo_p - (k_ass*AngII_Fl_Tb_exo*AT1R_memb_Tb ...
                                     - (k_diss + k_int)*AT1R_AngII_memb_Tb_exo);

    f(53) = AT1R_AngII_cell_Tb_exo_p - ((V_Tb_Fl/V_Tb_Pt_cell)*k_int*AT1R_AngII_memb_Tb_exo ...
                             + k_ass*AngII_cell_Tb_exo*AT1R_cell_Tb - k_diss*AT1R_AngII_cell_Tb_exo) ;

    f(54) = AngII_cell_Tb_exo_p - (k_diss*AT1R_AngII_cell_Tb_exo - k_ass*AngII_cell_Tb_exo*AT1R_cell_Tb ...
                                - k_lys*AngII_cell_Tb_exo);

    % Renal (blood) vasculature
    f(55) = AngII_Pv_exo_p - (((phi_RPF - phi_GFR - phi_L)/V_Pv)*AngII_circ_exo ...
                           + (phi_Pv/V_Pv)*AngII_Isf_Pt_exo...
                           - (((phi_RPF-2*phi_L-phi_U)/V_Pv)+v_II)*AngII_Pv_exo...
                           - k_ass*AT1R_memb_Pv*AngII_Pv_exo + k_diss*AT1R_AngII_memb_Pv_exo); 
               
    f(56) = AT1R_AngII_memb_Pv_exo_p - (k_ass*AT1R_memb_Pv*AngII_Pv_exo - k_diss*AT1R_AngII_memb_Pv_exo);
    
    % Total
    f(57) = AngII_T_exo - (V_Gl_Isf*AngII_Isf_Gl_exo + V_Pt_Isf*AngII_Isf_Pt_exo + V_Tb_Fl*AngII_Fl_Tb_exo...
                + V_Gl_Isf*AT1R_AngII_memb_Gl_exo + V_Pt_Isf*AT1R_AngII_memb_Pt_exo + V_Tb_Fl*AT1R_AngII_memb_Tb_exo...
                + V_Tb_Pt_cell*(AT1R_AngII_cell_Tb_exo + AT1R_AngII_cell_Pt_exo + AngII_cell_Tb_exo + AngII_cell_Pt_exo)...
                + V_Gl_cell*(AT1R_AngII_cell_Gl_exo + AngII_cell_Gl_exo) + V_Pv*(AngII_Pv_exo+AT1R_AngII_memb_Pv_exo));

    %%% --- Total
    
    % Systemic
    f(58) = AngII_circ_tot - (AngII_circ + AngII_circ_exo);
    f(59) = AT1R_AngII_memb_circ_tot - (AT1R_AngII_memb_circ + AT1R_AngII_memb_circ_exo);
    
    % Glomerular
    f(60) = AngII_Isf_Gl_tot - (AngII_Isf_Gl + AngII_Isf_Gl_exo);
    f(61) = AT1R_AngII_memb_Gl_tot - (AT1R_AngII_memb_Gl + AT1R_AngII_memb_Gl_exo);
    f(62) = AT1R_AngII_cell_Gl_tot - (AT1R_AngII_cell_Gl + AT1R_AngII_cell_Gl_exo);
    f(63) = AngII_cell_Gl_tot - (AngII_cell_Gl + AngII_cell_Gl_exo);
    
    % Peritubular
    f(64) = AngII_Isf_Pt_tot - (AngII_Isf_Pt + AngII_Isf_Pt_exo);
    f(65) = AT1R_AngII_memb_Pt_tot - (AT1R_AngII_memb_Pt + AT1R_AngII_memb_Pt_exo);
    f(66) = AT1R_AngII_cell_Pt_tot - (AT1R_AngII_cell_Pt + AT1R_AngII_cell_Pt_exo);
    f(67) = AngII_cell_Pt_tot - (AngII_cell_Pt + AngII_cell_Pt_exo);
    
    % Tubular
    f(68) = AngII_Fl_Tb_tot - (AngII_Fl_Tb + AngII_Fl_Tb_exo);
    f(69) = AT1R_AngII_memb_Tb_tot - (AT1R_AngII_memb_Tb + AT1R_AngII_memb_Tb_exo);
    f(70) = AT1R_AngII_cell_Tb_tot - (AT1R_AngII_cell_Tb + AT1R_AngII_cell_Tb_exo);
    f(71) = AngII_cell_Tb_tot - (AngII_cell_Tb + AngII_cell_Tb_exo);
    
    % Renal (blood) vasculature
    f(72) = AngII_Pv_tot - (AngII_Pv + AngII_Pv_exo);
    f(73) = AT1R_AngII_memb_Pv_tot - (AT1R_AngII_memb_Pv + AT1R_AngII_memb_Pv_exo);

    % Whole Kidney
    f(74) = AngII_T_tot - (AngII_T + AngII_T_exo);
end

end
