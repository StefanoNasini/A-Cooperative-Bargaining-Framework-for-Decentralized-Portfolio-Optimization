function [Grad] = Compute_CVaR_grad(z, R, Psi, alpha_B, d_CVaR)

%alpha_B=0.5;
B = exp( - (normal_quantile(alpha_B)^2)/2 ) /(alpha_B*sqrt(2*pi) );
n = length(R);
Grad = zeros(n,1);

for i = 1:n
    
    n_i = size(R{i},1);
    mu_R = mean(table2array(R{i}));
    mu_R_matrix = repmat(mu_R, n_i, 1); 
    c = norm((table2array(R{i}) - mu_R_matrix)*Psi{i});

    Grad(i) =  (  (1 - mu_R*Psi{i} + c )*B  )/( ( d_CVaR(i) - z(i)*(1 - mu_R*Psi{i} + c*B)*B  )   );
    
end


end