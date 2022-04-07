% ODE model for WinCF with antibiotic treatments

% parameters in the model:
%SetPara_WinCF_1D_Antibiotic; 

SetPara_WinCF_1D_Antibiotic_A_inhi_P; 


% create the matrix to hold values of all patients, all treatments, write
% to excel file at the end

data_PT = zeros(N_ratio, 2*N_Treat + 1); % each treat has two data point (theta_p, theta_f), and 1 for patient ID


% loop over the initial Psudomonas-fermenter biomass ratio
for I_ratio = 1: N_ratio
    
    %data_PT(I_ratio,1) = Sample_ID(I_ratio); % save the patiend_ID first
    
    data_PT(I_ratio,1) = lambda; % only for patient 12, save inhibition parameter lambda
    
    % set the initial biomass ratio of Pseudomonas to fermenters, assume it
    % is uniform on an interval of length 0.8
    theta_p_0 = Ratio_pf_ini(I_ratio,1)/0.8;
    theta_f_0 = Ratio_pf_ini(I_ratio,2)/0.8;
    
%     fname = sprintf('./data/PSpu_%d_pf_ratio_diff_treat',Sample_ID(I_ratio));
%     fid_50 = fopen(fname,'w');  % file to save the final biomass ratio by different treatments with same initial biomass ratio
%     
    for I_Treat = 1: N_Treat
        
        % Set the treatment combination
        Ind_Tp = Treat_ID(I_Treat,1);  % Indicator for antibiotic only kills Pseudomona
        Ind_Tf = Treat_ID(I_Treat,2);  % Indicator for antibiotic only kills anarobes
        Ind_Tw = Treat_ID(I_Treat,3);  % Indicator for wide spectrum antibiotic
        
        %initialize
        [theta_f_n, theta_p_n, F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, Tf_n, Tp_n, Tw_n] = Initialize_WinCF_1D_Antibiotic();
        
        fname = sprintf('./data/PSpu_%d_%s_theta_f',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid1 = fopen(fname,'w');  % file for fermentation bacteria population
        
        fname = sprintf('./data/PSpu_%d_%s_theta_p',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid2 = fopen(fname,'w'); % file for Pseudomonas aeruginosa population
        
        fname = sprintf('./data/PSpu_%d_%s_F',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid3 = fopen(fname,'w');      % file for chemical F
        
        fname = sprintf('./data/PSpu_%d_%s_P',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid4 = fopen(fname,'w');   % file for chemical P
        
        fname = sprintf('./data/PSpu_%d_%s_I',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid5 = fopen(fname,'w');      % file for I
        
        fname = sprintf('./data/PSpu_%d_%s_SO',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid6 = fopen(fname,'w');   % file for Oxygen
        
        fname = sprintf('./data/PSpu_%d_%s_SN',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid7 = fopen(fname,'w');   % file for Nitrate
        
        fname = sprintf('./data/PSpu_%d_%s_SA',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid8 = fopen(fname,'w');      % file for Amino-acid
        
        fname = sprintf('./data/PSpu_%d_%s_SG',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid9 = fopen(fname,'w');   % file for sugar
        
        fname = sprintf('./data/PSpu_%d_%s_Tf',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid10 = fopen(fname,'w');   % file for antibiotic only kills anarobes
        
        fname = sprintf('./data/PSpu_%d_%s_Tp',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid11 = fopen(fname,'w');   % file for antibiotic only kills Pseudomonas
        
        fname = sprintf('./data/PSpu_%d_%s_Tw',Sample_ID(I_ratio),Treat_Names(I_Treat));
        fid12 = fopen(fname,'w');   % file for wide spectrum antibiotic
        
        fprintf(fid1, '%12.10e ', theta_f_n);
        fprintf(fid1, '\n');
        
        fprintf(fid2, '%12.10e ', theta_p_n);
        fprintf(fid2, '\n');
        
        fprintf(fid3, '%12.10e ', F_n);
        fprintf(fid3, '\n');
        
        fprintf(fid4, '%12.10e ', P_n);
        fprintf(fid4, '\n');
        
        fprintf(fid5, '%12.10e ', I_n);
        fprintf(fid5, '\n');
        
        fprintf(fid6, '%12.10e ', SO_n);
        fprintf(fid6, '\n');
        
        fprintf(fid7, '%12.10e ', SN_n);
        fprintf(fid7, '\n');
        
        fprintf(fid8, '%12.10e ', SA_n);
        fprintf(fid8, '\n');
        
        fprintf(fid9, '%12.10e ', SG_n);
        fprintf(fid9, '\n');
        
        fprintf(fid10, '%12.10e ', Tf_n);
        fprintf(fid10, '\n');
        
        fprintf(fid11, '%12.10e ', Tp_n);
        fprintf(fid11, '\n');
        
        fprintf(fid12, '%12.10e ', Tw_n);
        fprintf(fid12, '\n');
        
        for nt = 1:NT  % loop time
            
            % Solve F, P, I, SO, SN, SA, SG, Tf, Tp, Tw at t_(n + 1), interpolate their value at t_(n + 1/)
            % Use  RK method for calculating theta_f and theta_p at t_(n+1)
            
            % F, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 1);
            
            [d,e,f,rhs] = Build_RD_tridiag(F_n, D_F, r_imp, r_exp, BC_flag_F, 0, 0, dt);
            
            F_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % P, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 2);
            
            [d,e,f,rhs] = Build_RD_tridiag(P_n, D_P, r_imp, r_exp, BC_flag_P, 0, 0, dt);
            
            P_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % I, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 3);
            
            [d,e,f,rhs] = Build_RD_tridiag(I_n, D_I, r_imp, r_exp, BC_flag_I, 0, 0, dt);
            
            I_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % SO, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 4);
            
            [d,e,f,rhs] = Build_RD_tridiag(SO_n, D_O, r_imp, r_exp, BC_flag_SO, 0, SO_BC1, dt);
            
            dSO_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % SN, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 5);
            
            [d,e,f,rhs] = Build_RD_tridiag(SN_n, D_N, r_imp, r_exp, BC_flag_SN, 0, 0, dt);
            
            SN_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % SA , calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 6);
            
            [d,e,f,rhs] = Build_RD_tridiag(SA_n, D_A, r_imp, r_exp, BC_flag_SA, 0, 0, dt);
            
            SA_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % SG, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 7);
            
            [d,e,f,rhs] = Build_RD_tridiag(SG_n, D_G, r_imp, r_exp, BC_flag_SG, 0, 0, dt);
            
            SG_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % Tf, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 8);
            
            [d,e,f,rhs] = Build_RD_tridiag(Tf_n, D_Tf, r_imp, r_exp, BC_flag_Tf, 0, 0, dt);
            
            Tf_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % Tp, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 9);
            
            [d,e,f,rhs] = Build_RD_tridiag(Tp_n, D_Tp, r_imp, r_exp, BC_flag_Tp, 0, 0, dt);
            
            Tp_np1 = Tridiag_Solver(d,e,f, rhs);
            
            % Tw, calculate the implicit and explicit reaction terms first
            [r_imp, r_exp] = Build_r_imp_exp(F_n, P_n, I_n, SO_n, SN_n, SA_n, SG_n, theta_f_n, theta_p_n, 10);
            
            [d,e,f,rhs] = Build_RD_tridiag(Tw_n, D_Tw, r_imp, r_exp, BC_flag_Tw, 0, 0, dt);
            
            Tw_np1 = Tridiag_Solver(d,e,f, rhs);
            
            
            % Impose Dirichlet BC on oxygen
            SO_np1 = SO_n;
            SO_np1(1:Nx) = dSO_np1;
            
            
            F_nph = .5*(F_n + F_np1);
            P_nph = .5*(P_n + P_np1);
            I_nph = .5*(I_n + I_np1);
            SO_nph = .5*(SO_n + SO_np1);
            SN_nph = .5*(SN_n + SN_np1);
            SA_nph = .5*(SA_n + SA_np1);
            SG_nph = .5*(SG_n + SG_np1);
            Tf_nph = .5*(Tf_n + Tf_np1);
            Tp_nph = .5*(Tp_n + Tp_np1);
            Tw_nph = .5*(Tw_n + Tw_np1);
            
            
            % use Runge Kutta method to solve the ODEs for theta_f and theta_p
            
            U = [theta_f_n, theta_p_n];
            
            Y1 = U;
            
            fY1 = FP_RHS_WinCF_1D_Antibiotic(Y1,F_n,P_n,I_n,SO_n,SN_n,SA_n,SG_n,Tf_n,Tp_n,Tw_n);
            Y2 = U + .5*dt*fY1;
            
            fY2 = FP_RHS_WinCF_1D_Antibiotic(Y2,F_nph,P_nph,I_nph,SO_nph,SN_nph,SA_nph,SG_nph,Tf_nph,Tp_nph,Tw_nph);
            Y3 = U + .5*dt*fY2;
            
            fY3 = FP_RHS_WinCF_1D_Antibiotic(Y3,F_nph,P_nph,I_nph,SO_nph,SN_nph,SA_nph,SG_nph,Tf_nph,Tp_nph,Tw_nph);
            Y4 = U +  dt*fY3;
            
            fY4 = FP_RHS_WinCF_1D_Antibiotic(Y4,F_np1,P_np1,I_np1,SO_np1,SN_np1,SA_np1,SG_np1,Tf_np1,Tp_np1,Tw_np1);
            %update U
            U = U + (dt/6)*(fY1+ 2.0*fY2 + 2.0*fY3 + fY4);
            
            theta_f_n = U(:,1);
            theta_p_n = U(:,2);
            
            % update the chemicals
            F_n = F_np1;
            P_n = P_np1;
            I_n = I_np1;
            SO_n =  SO_np1;
            SN_n = SN_np1;
            SA_n = SA_np1;
            SG_n = SG_np1;
            Tf_n = Tf_np1;
            Tp_n = Tp_np1;
            Tw_n = Tw_np1;
            
            
            if(Save_Inter_flag == 1 && mod(nt,FS_interval) == 0)
                
                fprintf(fid1, '%12.10e ', theta_f_n);
                fprintf(fid1, '\n');
                
                fprintf(fid2, '%12.10e ', theta_p_n);
                fprintf(fid2, '\n');
                
                fprintf(fid3, '%12.10e ', F_n);
                fprintf(fid3, '\n');
                
                fprintf(fid4, '%12.10e ', P_n);
                fprintf(fid4, '\n');
                
                fprintf(fid5, '%12.10e ', I_n);
                fprintf(fid5, '\n');
                
                fprintf(fid6, '%12.10e ', SO_n);
                fprintf(fid6, '\n');
                
                fprintf(fid7, '%12.10e ', SN_n);
                fprintf(fid7, '\n');
                
                fprintf(fid8, '%12.10e ', SA_n);
                fprintf(fid8, '\n');
                
                fprintf(fid9, '%12.10e ', SG_n);
                fprintf(fid9, '\n');
                
                fprintf(fid10, '%12.10e ', Tf_n);
                fprintf(fid10, '\n');
                
                fprintf(fid11, '%12.10e ', Tp_n);
                fprintf(fid11, '\n');
                
                fprintf(fid12, '%12.10e ', Tw_n);
                fprintf(fid12, '\n');
                
                
            end
            
        end
        
        if(Save_Inter_flag == 0) % save the final results if not saving the intermediate steps
            
            fprintf(fid1, '%12.10e ', theta_f_n);
            fprintf(fid1, '\n');
            
            fprintf(fid2, '%12.10e ', theta_p_n);
            fprintf(fid2, '\n');
            
            fprintf(fid3, '%12.10e ', F_n);
            fprintf(fid3, '\n');
            
            fprintf(fid4, '%12.10e ', P_n);
            fprintf(fid4, '\n');
            
            fprintf(fid5, '%12.10e ', I_n);
            fprintf(fid5, '\n');
            
            fprintf(fid6, '%12.10e ', SO_n);
            fprintf(fid6, '\n');
            
            fprintf(fid7, '%12.10e ', SN_n);
            fprintf(fid7, '\n');
            
            fprintf(fid8, '%12.10e ', SA_n);
            fprintf(fid8, '\n');
            
            fprintf(fid9, '%12.10e ', SG_n);
            fprintf(fid9, '\n');
            
            fprintf(fid10, '%12.10e ', Tf_n);
            fprintf(fid10, '\n');
            
            fprintf(fid11, '%12.10e ', Tp_n);
            fprintf(fid11, '\n');
            
            fprintf(fid12, '%12.10e ', Tw_n);
            fprintf(fid12, '\n');
        end
        
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
        fclose(fid4);
        fclose(fid5);
        fclose(fid6);
        fclose(fid7);
        fclose(fid8);
        fclose(fid9);
        fclose(fid10);
        fclose(fid11);
        fclose(fid12);
        % calculate the final total Pseudomonas and fermenter biomass, and their percentages, and save it
        
        total_p = sum(theta_p_n)*0.8/100;
        total_f = sum(theta_f_n)*0.8/100;
        perct_p = total_p/(total_p + total_f);
        perct_f = total_f/(total_p + total_f);
%         fprintf(fid_50,' %12.10e  %12.10e     %12.10e  %12.10e ', total_p, total_f, perct_p, perct_f);
%         fprintf(fid_50, '\n');
        data_PT(I_ratio,2*I_Treat) = total_p;
        data_PT(I_ratio,2*I_Treat+1) = total_f;
        
    end    
end


%fname = 'theta_pf_diff_ratio_treat.xlsx';

%only for patient 12, two treatment NT and Tf, with different lambda values
%(different inhibition strength)
fname = sprintf('theta_pf_P12_lambda_%1.2f.xlsx',lambda);


writematrix(data_PT,fname);








