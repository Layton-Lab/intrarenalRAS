function f = compute_regional_concentrations(sim)

% Volumes
V_Gl_Isf = 0.0019;
V_Pt_Isf = 0.0236;
V_Gl_cell = 0.0019;
V_Pv = 0.085;
V_Tb_Fl = 0.102;
V_Tb_Pt_cell = 0.294;

% Indices
Gl_Ce1_endo_i = 12; Gl_Ce1_exo_i = 46;
Gl_Cell_endo_i = 13; Gl_Cell_exo_i = 47;
Gl_Cs1_endo_i = 11; Gl_Cs1_exo_i = 45;
Gl_Isf_endo_i = 9; Gl_Isf_exo_i = 44;

Pt_Ce1_endo_i = 19; Pt_Ce1_exo_i = 50;
Pt_Cell_endo_i = 20; Pt_Cell_exo_i = 51;
Pt_Cs1_endo_i = 18; Pt_Cs1_exo_i = 49;
Pt_Isf_endo_i = 17; Pt_Isf_exo_i = 48;

Tb_Ce1_endo_i = 26; Tb_Ce1_exo_i = 54;
Tb_Cell_endo_i = 27; Tb_Cell_exo_i = 55;
Tb_Cs1_endo_i = 25; Tb_Cs1_exo_i = 53;
Tb_Fl_endo_i = 24; Tb_Fl_exo_i = 52;

Pv_Cs1_endo_i = 32; Pv_Cs1_exo_i = 57;
Pv_endo_i = 31; Pv_exo_i = 56;

T_endo_i = 35; T_exo_i = 58;

% Total intrarenal [Ang II]
T = sim(T_endo_i,:) + sim(T_exo_i,:);

Gl_intra_endo = V_Gl_cell*(sim(Gl_Ce1_endo_i,:) + sim(Gl_Cell_endo_i,:));
Gl_memb_endo = V_Gl_Isf*sim(Gl_Cs1_endo_i,:);
Gl_Isf_endo = V_Gl_Isf*sim(Gl_Isf_endo_i,:);
Pt_intra_endo = V_Tb_Pt_cell*(sim(Pt_Ce1_endo_i,:) + sim(Pt_Cell_endo_i,:));
Pt_memb_endo = V_Pt_Isf*sim(Pt_Cs1_endo_i,:);
Pt_Isf_endo = V_Pt_Isf*sim(Pt_Isf_endo_i,:);
Tb_intra_endo = V_Tb_Pt_cell*(sim(Tb_Ce1_endo_i,:) + sim(Tb_Cell_endo_i,:));
Tb_memb_endo = V_Tb_Fl*sim(Tb_Cs1_endo_i,:);
Tb_Fl_endo = V_Tb_Fl*sim(Tb_Fl_endo_i,:);
Pv_memb_endo = V_Pv*sim(Pv_Cs1_endo_i,:);
Pv_endo = V_Pv*sim(Pv_endo_i,:);

Gl_intra_exo = V_Gl_cell*(sim(Gl_Ce1_exo_i,:) + sim(Gl_Cell_exo_i,:));
Gl_memb_exo = V_Gl_Isf*sim(Gl_Cs1_exo_i,:);
Gl_Isf_exo = V_Gl_Isf*sim(Gl_Isf_exo_i,:);
Pt_intra_exo = V_Tb_Pt_cell*(sim(Pt_Ce1_exo_i,:) + sim(Pt_Cell_exo_i,:));
Pt_memb_exo = V_Pt_Isf*sim(Pt_Cs1_exo_i,:);
Pt_Isf_exo = V_Pt_Isf*sim(Pt_Isf_exo_i,:);
Tb_intra_exo = V_Tb_Pt_cell*(sim(Tb_Ce1_exo_i,:) + sim(Tb_Cell_exo_i,:));
Tb_memb_exo = V_Tb_Fl*sim(Tb_Cs1_exo_i,:);
Tb_Fl_exo = V_Tb_Fl*sim(Tb_Fl_exo_i,:);
Pv_memb_exo = V_Pv*sim(Pv_Cs1_exo_i,:);
Pv_exo = V_Pv*sim(Pv_exo_i,:)./T;

intra_endo = Gl_intra_endo + Pt_intra_endo + Tb_intra_endo;
memb_endo = Gl_memb_endo + Pt_memb_endo + Tb_memb_endo + Pv_memb_endo;
isf_endo = Gl_Isf_endo + Pt_Isf_endo + Tb_Fl_endo + Pv_endo;

intra_exo = Gl_intra_exo + Pt_intra_exo + Tb_intra_exo;
memb_exo = Gl_memb_exo + Pt_memb_exo + Tb_memb_exo + Pv_memb_exo;
isf_exo = Gl_Isf_exo + Pt_Isf_exo + Tb_Fl_exo + Pv_exo;

intra_tot = intra_exo + intra_endo;
memb_tot = memb_exo + memb_endo;
isf_tot = isf_exo + isf_endo;

endo_tot = intra_endo + memb_endo + isf_endo;
exo_tot = intra_exo + memb_exo + isf_exo;
tot = endo_tot + exo_tot;

f = [intra_endo; intra_exo; intra_tot;...
    memb_endo; memb_exo; memb_tot;...
    isf_endo; isf_exo; isf_tot;...
    endo_tot; exo_tot; tot];

end