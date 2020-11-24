function [OF] = Compute_s_moment_fun(z, R, Psi, s, d_moments)

n = length(R);

OF = 0;

for i = 1:n
    
    OF = OF - log( d_moments(i) - (z(i)^s)*sum((table2array(R{i})*Psi{i}).^s));
    
end

OF = real(OF);

end
 