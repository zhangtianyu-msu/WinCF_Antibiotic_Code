% Function to build the system of linear equations from discretizing the
% reaction-diffusion equation for all chemicals, the matrix is tridiagonal 

% Input:
%        C_n -> chemical concentration (vector of length Nx+1) at current time
%             step
%        D  -> diffusion coefficient
%        r_imp -> implicit reaction term, vector of length Nx + 1
%        r_exp -> explicit reaction term, vector of length Nx + 1
%        BC_flag -> Boundary condition flag, 1: No_flux at x = 0 and L; 
%                                            2: No_flux at x = 0, Dirichlet at x = L; 
%                                            3: Dirichlet at x = 0,No_flux at x = L;
%                                            4: Dirichlet at x = 0 and x = L
%        DBC0 -> the value of Dirichlet BC at x = 0 
%        DBC1 -> the value of Dirichlet BC at x = L, 
%        dt1 -> amount advanced in time axis, either dt or dt/2


% Output: d -> the main diagonal of the matrix (vector of length:  Nx + 1 for BC_flag = 1; length Nx for BC_flag = 2, 3; length Nx - 1 for BC_flag = 4)
%         e -> the subdiagonal of the matrix (vector of length Nx for BC_flag = 1; length Nx - 1 for BC_flag = 2, 3; length Nx - 2 for BC_flag = 4)
%         f -> the superdiagonal of the matrix (vector of length Nx for BC_flag = 1; length Nx - 1 for BC_flag = 2, 3; length Nx - 2 for BC_flag = 4)
%         rhs -> the right hand side of the linear system (vector of length:  Nx + 1 for BC_flag = 1; length Nx for BC_flag = 2, 3; length Nx - 1 for BC_flag = 4)


function [d,e,f,rhs] = Build_RD_tridiag(C_n, D, r_imp, r_exp, BC_flag, DBC0, DBC1, dt1)

global  dx Nx 

eta = dt1*D/(dx*dx);

if BC_flag == 1 % No_flux at x = 0 and L;
        
    n = Nx + 1;  % size of the tridiagonal linear system
    
    d = (1 + 2*eta)*ones(n, 1) - dt1*r_imp;
    e1 = -eta*ones(n, 1);
    f1 = -eta*ones(n, 1);
        
    rhs = C_n  + dt1*r_exp;
    
    f1(1) = f1(1) + e1(1); % no-flux BC at x = 0
    
    e1(n) = f1(n) + e1(n); % no-flux BC at x = L  
    
elseif BC_flag == 2 %  no flxu at x = 0, Dirichlet at x = L
    
    n = Nx;  % size of the tridiagonal linear system
    
    d = (1 + 2*eta)*ones(n, 1) - dt1*r_imp(1:Nx);
    e1 = -eta*ones(n, 1);
    f1 = -eta*ones(n, 1);
    
    % d = d + (dt1*r_a/Y_pO)*theta_p_n(1:Nx)./(K_O + SO_n(1:Nx));
    
    rhs = C_n(1:Nx) + dt1*r_exp(1:Nx);
    
    f1(1) = f1(1) + e1(1); % no-flux BC at x = 0
    
    rhs(n) = rhs(n) - f1(n)*DBC1; % Dirichlet at x = L
    
elseif BC_flag == 3 % Dirichlet at x = 0,No_flux at x = L;
    
    n = Nx;  % size of the tridiagonal linear system
    
    d = (1 + 2*eta)*ones(n, 1) - dt1*r_imp(2:Nx+1);
    e1 = -eta*ones(n, 1);
    f1 = -eta*ones(n, 1);
    
    % d = d + (dt1*r_a/Y_pO)*theta_p_n(1:Nx)./(K_O + SO_n(1:Nx));
    
    rhs = C_n(2:Nx+1) + dt1*r_exp(2:Nx+1);
    
    rhs(1) = rhs(1) - e1(1)*DBC0; % Dirichlet at x = 0
    
    e1(n) = f1(n) + e1(n); % no-flux BC at x = L
    
    
elseif BC_flag == 4 % Dirichlet at x = 0 and x = L
    
    n = Nx - 1;  % size of the tridiagonal linear system
    
    d = (1 + 2*eta)*ones(n, 1) - dt1*r_imp(2:Nx);
    e1 = -eta*ones(n, 1);
    f1 = -eta*ones(n, 1);
    
    % d = d + (dt1*r_a/Y_pO)*theta_p_n(2:Nx)./(K_O + SO_n(2:Nx));
    
    rhs = C_n(2:Nx) + dt1*r_exp(2:Nx);
    
    rhs(1) = rhs(1) - e1(1)*DBC0; % Dirichlet at x = 0
    
    rhs(n) = rhs(n) - f1(n)*DBC1; % Dirichlet at x = L
    
    
    
else
    
    disp(' The BC_flag has invalid value in Build_SO, must be 1 or 2! ');
    
end

        

e = e1(2:n);

f = f1(1:n-1);

end



