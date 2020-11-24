function [OF] = Compute_CVaR_fun(z, R, Psi, alpha_B, d_CVaR)

%R=follower;
%z=x;

%alpha_B=0.5;
B = exp( - (normal_quantile(alpha_B)^2)/2 ) /(alpha_B*sqrt(2*pi) );
n = length(R);
OF = 0;


for i = 1:n

    n_i = size(R{i},1);
    mu_R = mean(table2array(R{i}));
    mu_R_matrix = repmat(mu_R, n_i, 1); 
    c = norm((table2array(R{i}) - mu_R_matrix)*Psi{i});
    
    OF = OF - log( d_CVaR(i) - z(i)*(1 + mu_R*Psi{i} - c*B));
             
end

OF = real(OF);

end
 