%function to initialize the unknowns
function [theta_f_n, theta_p_n, F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, Tf_n, Tp_n, Tw_n] = Initialize_WinCF_1D_Antibiotic()

global Nx  theta_f_0 theta_p_0

global F_IC0 P_IC0 I_IC0 SO_BC1 SN_IC0 SA_IC0 SG_IC0 Tf_IC0 Tp_IC0 Tw_IC0

theta_f_n = theta_f_0*ones(Nx+1,1);

theta_p_n = theta_p_0*ones(Nx+1,1);

F_n = F_IC0*ones(Nx+1,1);

P_n = P_IC0*ones(Nx+1,1);

I_n = I_IC0*ones(Nx+1,1);

SO_n = SO_BC1*ones(Nx+1,1);

SN_n = SN_IC0*ones(Nx+1,1);

SA_n = SA_IC0*ones(Nx+1,1);

SG_n = SG_IC0*ones(Nx+1,1);

Tf_n = Tf_IC0*ones(Nx+1,1);

Tp_n = Tp_IC0*ones(Nx+1,1);

Tw_n = Tw_IC0*ones(Nx+1,1);

end

 
 