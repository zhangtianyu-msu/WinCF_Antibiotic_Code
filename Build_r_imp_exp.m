% function to build the implicit and explicit part of the reaction terms
% for the reaction-diffusion equations
% C_t = D*C_xx + C*r_imp + r_exp
% Input: F, P, I, SO, SN, SA, SG, theta_f, theta_p  ----
%        currently value of all chemical concentrations and bacteria
%        density , all are vectors of length Nx + 1
%        flag: 1 = F, 2 = P, 3 = I, 4 = SO, 5 = SN, 6 = SA, 7 = SG, 8 = Tf, 9 = Tp, 10 = Tw 
% Output: r_imp -> implicit reaction term
%         r_exp -> explicity reatciotn term

function [r_imp, r_exp] = Build_r_imp_exp(F, P, I, SO, SN, SA, SG, theta_f, theta_p, flag)

global K_f K_p K_O K_N K_G K_A  Y_pPo Y_pPn  Y_pO Y_pN Y_pIo Y_pIn Y_pA  Y_fF Y_fG d_O d_p %d_F
global d_Tp d_Tf d_pw d_fw Y_pT Y_fT Y_pw Y_fw % dimensionless parameters, global variables

global mu_f mu_pa mu_pn  %growth rates of bacteria

global Nx

global beta0 beta1 beta2 C_0 Ind_Tp Ind_Tf Ind_Tw  % control parameters 
 
g0 = @(x) mu_f*(1-(2/pi)*atan(beta0*x));
g1 = @(x) K_f*(1 - 0.9*(2/pi)*atan(beta1*x));
g2 = @(x) K_p*(0.9*(0.5 - (1/pi)*atan(beta2*x)) + 0.1);

r_imp = zeros(Nx+1,1);
r_exp = zeros(Nx+1,1);

nonneg_g1_I = 1; %max(0,1 - theta_f./g1(I)); % non-negative logistic term associated with I, used in production of chemical to avoid negative production rate
nonneg_g2_pH = 1; %max(0,1 - theta_p./g2(F - P + C_0)); % non-negative logistic term associated with pH, used in production of chemical to avoid negative production rate

if flag == 1 % F
    
    % no implicit reaction term
    %r_imp = -d_F*ones(Nx+1,1);  % This is for the decay of F
    
    r_exp = (1/Y_fF)*g0(SO).*(SG./(K_G + SG)).*theta_f.*(nonneg_g1_I);
    
%    keyboard;
    
      
elseif flag == 2 % P
    
    % no implicit reaction term
    
    r_exp = ( (mu_pa/Y_pPo)*SO./(K_O + SO) + (mu_pn*exp(-d_p*SO)/Y_pPn).*SN./(K_N + SN)).*(SA./(K_A + SA)).*theta_p.*(nonneg_g2_pH);
    
elseif flag == 3 % I 
    
    % no implicit reaction term
    
    r_exp = ( (mu_pa/Y_pIo)*SO./(K_O + SO) + (mu_pn*exp(-d_p*SO)/Y_pIn).*SN./(K_N + SN)).*(SA./(K_A + SA)).*theta_p.*(nonneg_g2_pH);
    
elseif flag == 4 % SO
    
    % no explicit reaction term
    
    r_imp = -theta_p.*( ((mu_pa/Y_pO)./(K_O + SO)).*(SA./(K_A + SA)).*(nonneg_g2_pH) + (d_O/Y_pO)./(K_O + SO) );
    
    
elseif flag == 5 % SN
    
    % no explicit reaction term
    
    r_imp = - (mu_pn*exp(-d_p*SO)/Y_pN)./(K_N + SN).*(SA./(K_A + SA)).*theta_p.*(nonneg_g2_pH);
    
    
elseif flag == 6 % SA
    
    % no explicit reaction term
    
    r_imp = -(1/Y_pA)*(mu_pa*SO./(K_O + SO) + mu_pn*exp(-d_p*SO).*SN./(K_N + SN))./(K_A + SA).*theta_p.*(nonneg_g2_pH);
    
elseif flag == 7 % SG
    
    % no explicit reaction term
    
    r_imp = -(1/Y_fG)*g0(SO)./(K_G + SG).*theta_f.*(nonneg_g1_I);
    
    
elseif flag == 8 % Tf
    
    % no explicit reaction term
    
    r_imp = - (Ind_Tf/Y_fT)*d_Tf*theta_f;   

elseif flag == 9 % Tp
    
    % no explicit reaction term
    
    r_imp = - (Ind_Tp/Y_pT)*d_Tp*theta_p;   
      
elseif flag == 10 % Tw
    
    % no explicit reaction term
    
    r_imp = - Ind_Tw*(d_fw*theta_f/Y_fw + d_pw*theta_p/Y_pw); 
    
  
else
   
     fprintf(1,'\n flag must be between 1 and 10, error in Build_r_imp_exp.m !!! \n');


end

