% function to define the RHS of the ODEs 
% U' = F(U,F,P,I,SO,SN,SA,SG,Tf,Tp,Tw), here U = (theta_f; theta_p) 
% array of size (Nx+1) by 2


function Y = FP_RHS_WinCF_1D_Antibiotic(U,F,P,I,SO,SN,SA,SG,Tf,Tp,Tw)

global mu_f mu_pa mu_pn K_f K_p K_O K_N  K_A K_G  beta0 beta1 beta2 C_0 d_p 
global d_Tp d_Tf d_pw d_fw Ind_Tp Ind_Tf Ind_Tw 
global lambda % parameter to adjust pH inhibition on Pseudomonas



g0 = @(x) mu_f*(1-(2/pi)*atan(beta0*x));
%g1 = @(x) K_f*(1 - 0.9*(2/pi)*atan(beta1*x));
g1 = @(x) K_f*(1 - 0.95*(2/pi)*atan(beta1*x));
%g2 = @(x) K_p*(0.9*(0.5 - (1/pi)*atan(beta2*x)) + 0.1);
%g2 = @(x) K_p*(0.95*(0.5 - (1/pi)*atan(beta2*x)) + 0.05);
g2 = @(x) K_p*((1 - lambda)*(0.5 - (1/pi)*atan(beta2*x)) + lambda);

Y = zeros(size(U));

tmp1 = g0(SO);
tmp2 = (SG./(K_G + SG)).*U(:,1).*(1-U(:,1)./g1(I));
tmp1(tmp2 < 0) = 1; % only use the oxygen inhibition when the growth rate of fermenters is positive

%Y(:,1) = g0(SO).*(SG./(K_G + SG)).*U(:,1).*(1-U(:,1)./g1(I)) - (Ind_Tf*d_Tf*Tf + Ind_Tw*d_fw*Tw).*U(:,1);  
Y(:,1) = tmp1.*tmp2 - (Ind_Tf*d_Tf*Tf + Ind_Tw*d_fw*Tw).*U(:,1);  

Y(:,2) = (mu_pa*(SO./(K_O + SO)) + exp(-d_p*SO)*mu_pn.*(SN./(K_N + SN))).*(SA./(K_A + SA)).*U(:,2).*(1 - U(:,2)./g2(F - P + C_0)) - (Ind_Tp*d_Tp*Tp + Ind_Tw*d_pw*Tw).*U(:,2);



end





