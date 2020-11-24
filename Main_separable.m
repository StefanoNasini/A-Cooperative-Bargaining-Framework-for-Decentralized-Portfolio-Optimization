
%--------------------------------------------------------------------------
% Install the non-linear convex solver
% cd 'C:\Users\s.nasini\Dropbox\Francisco+Stefano+Rabia\Programs\cvx'
% cvx_setup
%--------------------------------------------------------------------------

fclose('all');
close all;
clear all;
clc;

delete('.\output\*.*');

addpath('C:\Users\s.nasini\Dropbox\Francisco+Stefano+Rabia\1. Data\ProcessedData') 
addpath('C:\Users\s.nasini\Dropbox\cvx')

addpath('C:\Users\franc\Dropbox\Decentralized_Investments_Multiple_Intermediaries\1. Data\FollowersTestData\Old_followers_300assets') 
addpath('C:\Users\franc\Dropbox\Decentralized_Investments_Multiple_Intermediaries\1. Data\ProcessedData') 

nn = 73;

payoff_set = [0 1];                 % If 0 Moment, if 1 CVaR
alpha_set = [0.8 0.9];
J_set = [2 4];
rho_set = [0.05 0.25];               % Minimum return accepted

delta = 0.1;

instance = 1;

fid_summary = fopen('.\output\SummaryTable.log','a+');
fprintf(fid_summary,'Instance Year Payoff j/a rho gap step_k Leader_OF_opt Leader_OF_x0 Leader_OF_x1 top_opt top_x1 z_max z_min z_mean top_z Iterations CPU time\r\n');
            
for payoff = payoff_set

    if payoff == 0 
        JA_set = J_set;
    else
        JA_set = alpha_set;
    end
        
    for JA = JA_set

        for rho = rho_set
            
            data = cell(1,nn); 
            follower = cell(1,nn); 
            follower_mat = cell(1,nn);
            Year = 1999:2014;
            year_start = Year(1);
            year_end = Year(end);
            number_of_years = year_end-year_start + 1; 
            ActiveMarkets = cell(1,number_of_years); 
            
            Tollerance = 1e-5;                  %Tolerance between solution x and x+1

            for YY = 1:number_of_years

                ii = 0;
                for i = 1:nn
                    
                  FOLLOWER_FILE_NAME = ['Follower_' num2str(i,'%d') '.csv']; 
                  data{i} = readtable(FOLLOWER_FILE_NAME ); 
                  %follower{i} = data{i}(data{i}.Year == Year(YY), :);
                  [rown coln] = size(data{i});
                  follower{i} = table2cell(data{i}((12*(YY-1) + 1):(12*(YY-1) + 12), 4:coln));
                  
                  follower_mat{i} = [];
                  for cc = 1:12
                        follower_mat{i} = [follower_mat{i}; cell2mat(follower{i}(cc,:))];
                  end
                  
                  [I J] = find(follower_mat{i} <= -1000);
                  J = unique(J);
                  
                  follower_mat{i}(:,J) = [];
                  
                  if size(follower_mat{i},2) > 1      
                    ii = ii + 1;
                  end
                  
                  %[size(follower_mat{i},2) i ii]
                  %follower{i}.t = []; % Eliminate column t (not needed)
                  %follower{i}.Date = []; % Eliminate column Data (not needed)
                  %follower{i}.Year = []; % Eliminate column Year (not needed)
                                  
                end
                
                n = ii;
                R = cell(1,n); 
                                
                ii = 1;
                InactiveMarkets{YY} = [];
                for i = 1:nn
                    if size(follower_mat{i},2) > 1  
                        %[size(follower_mat{i},2) i ii]
                        R{ii} = array2table(follower_mat{i});
                        ActiveMarkets{YY} = [ActiveMarkets{YY} i];   
                        ii = ii + 1;
                    end
                end
                
                [OF, iterCount, tEnd, x, x0_initial_point, Psi, d] = Solve_separable(delta, JA, R, JA, rho, Tollerance, payoff);
                
                if payoff == 0 
                    OF_x0 = Compute_s_moment_fun(x0_initial_point, R, Psi, JA, d);
                
                    top_x1 = 1;
                    OF_x1_old = Inf;
                    for k = 1:n
                        x1 = zeros(1,n);
                        x1(k) = 1;
                        OF_x1 = min(OF_x1_old, Compute_s_moment_fun(x1, R, Psi, JA, d));
                        if OF_x1 < OF_x1_old
                            top_x1 = k;
                            OF_x1_old = OF_x1;
                        end
                    end
                        
                else
                    OF_x0 = Compute_CVaR_fun(x0_initial_point, R, Psi, JA, d);
                    
                    top_x1 = 1;
                    OF_x1_old = Inf;
                    for k = 1:n
                        x1 = zeros(1,n);
                        x1(k) = 1;
                        OF_x1 = min(OF_x1_old, Compute_CVaR_fun(x1, R, Psi, JA, d));
                        if OF_x1 < OF_x1_old
                            top_x1 = k;
                            OF_x1_old = OF_x1;
                        end
                    end
                end
    
                top_opt = find(x == max(x));
                
                XX{YY} = x';
                ValueOfx{YY}= Psi';
                ValueOfz{YY} = x';
                ValueOfz_initial{YY} = x0_initial_point'; 
                V{YY} = OF;
                
                % fid_table = fopen(['.\output\SummaryTable_instance_' num2str(instance,'%d') '_Year' num2str(Year(YY),'%d') '.log'],'a+');
                % [nrows, ncols] = size(x);
                % for i = 1:nrows
                %     fprintf(fid_table,[repmat(' %f\t ', 1, ncols) '\r\n'], round(x(i,:),3));
                % end
                % fprintf(fid_table,'\r\n');
                % fclose(fid_table); 

                z_max = find(x==max(x));
                
                fprintf(fid_summary,'%d %d %d %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %d %-8.2f %-8.2f \r\n', instance, Year(YY), payoff, JA, rho, Tollerance, OF, OF_x0, OF_x1, top_opt, top_x1, max(x), min(x), mean(x), z_max, iterCount, tEnd);

            end%year
      
            Make_plot;
            
            clear data follower Year year_start year_end number_of_years ;
            
            instance = instance + 1;
            
        end%rho
        
    end%JA
    
end%payoff

fclose(fid_summary);


