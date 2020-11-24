function [Psi] = SolveFollower_CVaR(R, alpha_B, rho)
  
n = length(R);
MeanReturns = cell(1,n);
%Cov = cell(1,n);
Psi = cell(1,n); 
%alpha_B=0.5;
B = exp( - (normal_quantile(alpha_B)^2)/2 ) /(alpha_B*sqrt(2*pi) );

%--------------------------------------------------------------------------
% Obtaining Psi_i
%--------------------------------------------------------------------------

for i=1:n
    
    n_i = size(R{i},1);
    MeanReturns{i} = repmat(mean(table2array(R{i})), n_i, 1); 
    n_k = size(R{i},2);
    xx0 = repmat(1/n_k,n_k,1);
    cvx_clear
    cvx_begin
    variable x(size(R{i},2),1); 
    c = (table2array(R{i}) - MeanReturns{i})*x;
    min( 1 - MeanReturns{i}*x  -  norm(c)*B );
    subject to% 
        sum(x) == 1;
        x'*MeanReturns{i}' >= rho*(xx0'*MeanReturns{i}');
        x >= 0;
    cvx_end
    x(isnan(x))=0;
    Psi{i} = x;
    
end

end