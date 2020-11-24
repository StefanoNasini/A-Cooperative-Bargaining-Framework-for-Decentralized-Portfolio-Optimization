
%--------------------------------------------------------------------------
% Aggregate plot for a given combination of (payoff, JA, rho)
%--------------------------------------------------------------------------
col_map = gray(nn);
%uncomment this lines to select different step size gray color
%c1_3 = [1:n]; 
%c1_3 = 3*c1_3;
%c1_3 = c1_3';
%AA = col_map(2:2:end,:);

for cc=1:nn
c1(cc) = strcat('[',  num2str(col_map(cc,1)),  {' '} ,num2str(col_map(cc,2)),  {' '} ,num2str(col_map(cc,3)),  ']' );
end 
c1=c1';

leader_OF = [];

z_percentages = zeros(nn,number_of_years);

for YY = 1:number_of_years
    z_percentages(ActiveMarkets{YY} ,YY) = XX{YY};
    leader_OF = [leader_OF V{YY}];
end

MyPlot = figure(1);
ax1 = axes('Parent', MyPlot );
hold(ax1, 'on');
b=bar(Year, z_percentages.','stacked');
for i = 1:length(b)
    b(i).FaceColor = c1{i};
end
p=plot(Year,(leader_OF-min(leader_OF))/max(1+(leader_OF-min(leader_OF))), 'bo--');
p.LineWidth = 2;
p.Color = 'red';
ylabel('Partition and Risk','FontSize',12,'FontWeight','bold','Color','r');
xlabel('Year','FontSize',12,'FontWeight','bold','Color','r')
ylim([0 1]);
xlim([1998 2015]);


hold(ax1, 'off');

if payoff == 0
    saveas(gcf,['.\output\Moments_instance_' num2str(instance,'%d') '.pdf']);
else
    saveas(gcf,['.\output\CVaR_instance_' num2str(instance,'%d') '.pdf']);
end

close all
close all hidden
close all force
