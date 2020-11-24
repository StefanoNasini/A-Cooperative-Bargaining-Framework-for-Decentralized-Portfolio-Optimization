function [Psi] = SolveFollower_s_moment(R, s, rho)
  
n = length(R);
MeanReturns = cell(1,n);
%Cov = cell(1,n);
Psi = cell(1,n); 

%--------------------------------------------------------------------------
% Obtaining Psi_i
%--------------------------------------------------------------------------

for i = 1:n
    
    MeanReturns{i} = mean(table2array(R{i})); 
    %Cov{i} = cov(table2array(R{i}));
    %[row, col] = find(isnan(MeanReturns{i}))
    n_k = size(R{i},2);
    xx0 = repmat(1/n_k,n_k,1);
    cvx_clear
    cvx_begin
    variable x(n_k,1); 
    %minimize(0.5*(x'*Cov{i}*x));  
    minimize( sum((table2array(R{i})*x).^s)) ;
    subject to% 
        sum(x) == 1;
        x'*MeanReturns{i}' >= rho*(xx0'*MeanReturns{i}');
        x >= 0;
    cvx_end
    x(isnan(x))=0;
    Psi{i} = x;
    
end

end


