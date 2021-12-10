function [comp_intra_tot, comp_memb_tot, comp_isf_tot, comp_tot] = get_distributions(sim,separate)

% Volumes
V_Gl_Isf = 0.0019;
V_Pt_Isf = 0.0236;
V_Gl_cell = 0.0019;
V_Pv = 0.147*0.58;
V_Tb_Fl = 0.102;
V_Tb_Pt_cell = 0.294;

% Indices
Gl_Ce1_endo_i = 12; Gl_Ce1_exo_i = 45;
Gl_Cell_endo_i = 13; Gl_Cell_exo_i = 46;
Gl_Cs1_endo_i = 11; Gl_Cs1_exo_i = 44;
Gl_Isf_endo_i = 9; Gl_Isf_exo_i = 43;

Pt_Ce1_endo_i = 19; Pt_Ce1_exo_i = 49;
Pt_Cell_endo_i = 20; Pt_Cell_exo_i = 50;
Pt_Cs1_endo_i = 18; Pt_Cs1_exo_i = 48;
Pt_Isf_endo_i = 17; Pt_Isf_exo_i = 47;

Tb_Ce1_endo_i = 26; Tb_Ce1_exo_i = 53;
Tb_Cell_endo_i = 27; Tb_Cell_exo_i = 54;
Tb_Cs1_endo_i = 25; Tb_Cs1_exo_i = 52;
Tb_Fl_endo_i = 24; Tb_Fl_exo_i = 51;

Pv_Cs1_endo_i = 32; Pv_Cs1_exo_i = 56;
Pv_endo_i = 31; Pv_exo_i = 55;

T_endo_i = 35; T_exo_i = 57;

if separate == true
    T = sim(T_endo_i,:) + sim(T_exo_i,:);
else 
    T = sim(T_endo_i,:);
end

Gl_intra_endo = V_Gl_cell*(sim(Gl_Ce1_endo_i,:) + sim(Gl_Cell_endo_i,:))./T;
Gl_memb_endo = V_Gl_Isf*sim(Gl_Cs1_endo_i,:)./T;
Gl_Isf_endo = V_Gl_Isf*sim(Gl_Isf_endo_i,:)./T;
Pt_intra_endo = V_Tb_Pt_cell*(sim(Pt_Ce1_endo_i,:) + sim(Pt_Cell_endo_i,:))./T;
Pt_memb_endo = V_Pt_Isf*sim(Pt_Cs1_endo_i,:)./T;
Pt_Isf_endo = V_Pt_Isf*sim(Pt_Isf_endo_i,:)./T;
Tb_intra_endo = V_Tb_Pt_cell*(sim(Tb_Ce1_endo_i,:) + sim(Tb_Cell_endo_i,:))./T;
Tb_memb_endo = V_Tb_Fl*sim(Tb_Cs1_endo_i,:)./T;
Tb_Fl_endo = V_Tb_Fl*sim(Tb_Fl_endo_i,:)./T;
Pv_memb_endo = V_Pv*sim(Pv_Cs1_endo_i,:)./T;
Pv_endo = V_Pv*sim(Pv_endo_i,:)./T;

if separate == true
    Gl_intra_exo = V_Gl_cell*(sim(Gl_Ce1_exo_i,:) + sim(Gl_Cell_exo_i,:))./T;
    Gl_memb_exo = V_Gl_Isf*sim(Gl_Cs1_exo_i,:)./T;
    Gl_Isf_exo = V_Gl_Isf*sim(Gl_Isf_exo_i,:)./T;
    Pt_intra_exo = V_Tb_Pt_cell*(sim(Pt_Ce1_exo_i,:) + sim(Pt_Cell_exo_i,:))./T;
    Pt_memb_exo = V_Pt_Isf*sim(Pt_Cs1_exo_i,:)./T;
    Pt_Isf_exo = V_Pt_Isf*sim(Pt_Isf_exo_i,:)./T;
    Tb_intra_exo = V_Tb_Pt_cell*(sim(Tb_Ce1_exo_i,:) + sim(Tb_Cell_exo_i,:))./T;
    Tb_memb_exo = V_Tb_Fl*sim(Tb_Cs1_exo_i,:)./T;
    Tb_Fl_exo = V_Tb_Fl*sim(Tb_Fl_exo_i,:)./T;
    Pv_memb_exo = V_Pv*sim(Pv_Cs1_exo_i,:)./T;
    Pv_exo = V_Pv*sim(Pv_exo_i,:)./T;

    Gl_intra = Gl_intra_endo + Gl_intra_exo;
    Gl_memb = Gl_memb_endo + Gl_memb_exo;
    Gl_Isf = Gl_Isf_endo + Gl_Isf_exo;
    Pt_intra = Pt_intra_endo + Pt_intra_exo;
    Pt_memb = Pt_memb_endo + Pt_memb_exo;
    Pt_Isf = Pt_Isf_endo + Pt_Isf_exo;
    Tb_intra = Tb_intra_endo + Tb_intra_exo;
    Tb_memb = Tb_memb_endo + Tb_memb_exo;
    Tb_Fl = Tb_Fl_endo + Tb_Fl_exo;
    Pv_memb = Pv_memb_endo + Pv_memb_exo;
    Pv = Pv_endo + Pv_exo;
else
    Gl_intra = Gl_intra_endo;
    Gl_memb = Gl_memb_endo;
    Gl_Isf = Gl_Isf_endo;
    Pt_intra = Pt_intra_endo;
    Pt_memb = Pt_memb_endo;
    Pt_Isf = Pt_Isf_endo;
    Tb_intra = Tb_intra_endo;
    Tb_memb = Tb_memb_endo;
    Tb_Fl = Tb_Fl_endo;
    Pv_memb = Pv_memb_endo;
    Pv = Pv_endo;  
end

Gl_tot = Gl_intra + Gl_memb + Gl_Isf;
Pt_tot = Pt_intra + Pt_memb + Pt_Isf;
Tb_tot = Tb_intra + Tb_memb + Tb_Fl;
Pv_tot = Pv_memb + Pv;

comp_intra_tot = [Gl_intra' Pt_intra' Tb_intra' zeros(size(Gl_intra'))];
comp_memb_tot = [Gl_memb' Pt_memb' Tb_memb' Pv_memb'];
comp_isf_tot = [Gl_Isf' Pt_Isf' Tb_Fl' Pv'];
comp_tot = [Gl_tot' Pt_tot' Tb_tot' Pv_tot'];

end