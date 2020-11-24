function [OF k tEnd x x0_initial_point Psi d] = Solve_separable(delta, alpha_B, R,  s, rho, Tollerance, payoff)

tStart = tic;    
n = length(R);

if payoff == 0
    
    Psi = SolveFollower_s_moment(R, s, rho);
    
    d_moments = ones(n,1);
    for i = 1:n
        %d_moments(i) = delta*sum((table2array(R{i})*Psi{i}).^s);
        %abs(d_moments(i) - 1/(size(R{i},2)^2))/max(abs(d_moments(i)), 1/(size(R{i},2)^2));
        %[d_moments(i) (1/(size(R{i},2)^2))]
        d_moments(i) = delta/(((s-1)^2)*size(R{i},2)^2);
    end

    d = d_moments;
    
else
    
    d_CVaR = ones(n,1);
    Psi = SolveFollower_CVaR(R, alpha_B, rho);
    
    B = exp( - (normal_quantile(alpha_B)^2)/2 ) /(alpha_B*sqrt(2*pi) );
    
    for i = 1:n
        
        n_i = size(R{i},1);
        mu_R = mean(table2array(R{i}));
        mu_R_matrix = repmat(mu_R, n_i, 1); 
        c = norm((table2array(R{i}) - mu_R_matrix)*Psi{i});
    
        %d_CVaR(i) = delta*(1 - mu_R*Psi{i} - c*B);
        d_CVaR(i) =  delta/(sqrt(size(R{i},2)));
        
    end
    
    d = d_CVaR;
    
end


%--------------------------------------------------------------------------
% Penality algorithm for the separable Knapsack Problem
%--------------------------------------------------------------------------

x = repmat(1/n,1,n);
x0 = x;
k = 0;
x_dif = 1;
x0_initial_point = x0;

while max(abs(x_dif)) > Tollerance 
    
        fprintf('----- Iteration %d --------------------------------------\r\n', k);
        
        if payoff == 0
            gradient = Compute_s_moment_grad(x, R, Psi, s, d_moments);
            gradient(isnan(gradient))= 0;
        else
            gradient = Compute_CVaR_grad(x, R, Psi, alpha_B, d_CVaR);
            gradient(isnan(gradient))= 0;
        end
        
        %beta = (20*n)/(k+1); 
        beta = (100*sqrt(n))/(k+1);
        gamma = max(1,sum(gradient.^2));
        alpha = beta/gamma;
        
        ff = exp(-alpha*(gradient));
        
        ff(isnan(gradient))= 0;
        
        fff = sum(x*ff);

        x = x.*(ff/fff)';

        
        
        x_dif = x0 - x;
        x0 = x;
        
	    k = k + 1;

end

tEnd = toc(tStart)  ; 


if payoff == 0
    OF = Compute_s_moment_fun(x, R, Psi, s, d_moments);
else
    OF = Compute_CVaR_fun(x, R, Psi, alpha_B, d_CVaR);
end
    
end