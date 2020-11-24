function [Grad] = Compute_s_moment_grad(z, R, Psi, s, d_moments)

n = length(R);

Grad = zeros(n,1);


for i = 1:n
    
    A_i = sum((table2array(R{i})*Psi{i}).^s);
    
    Grad(i) = -(s*(z(i)^(s-1))*A_i)/( d_moments(i) - (z(i)^s)*A_i);
    
end


end