%https://doi.org/10.1016/j.enconman.2008.09.003
%https://github.com/Raphael-Boichot/Tree-network-structure-generation-for-heat-conduction-by-cellular-automaton
%Raphaël BOICHOT initial VB code 2006, initial Matlab code 2009, published 2023
clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------Conditions for thermal science-----------------------------------
kp_k0=100;                      %kp/ko ratio, k0=1 by default
filling_ratio=0.3;              %filling ratio
heat_sink_temperature=298;      %self explanatory
delta_x=0.001;                  %step in x/y
p_vol=1e6;                      %surface power or volume density
starting_image='100x200.bmp';   %self explanatory
%---------Hyper parameters for the cellular automaton----------------------------
variation_rate=0.5;             %initial automaton variation rate, decreasing with epoch
verbose=0;                      %to suppress output other than final T max
min_ratio=0;                    %0% attraction by gradients only
max_ratio=0.15;                 %100% attraction by temperature only
steps=25;                       %searching steps
%--------------------------------------------------------------------------

ratios=linspace(min_ratio,max_ratio,steps);
a=zeros(1,length(ratios));
disp('Searching for best attraction parameter with small image')
parfor i=1:length(ratios)
    a(i)=automate_cell_direct(ratios(i),kp_k0,filling_ratio,heat_sink_temperature,delta_x,p_vol,variation_rate,starting_image,verbose);
end

% It's Gemini time!________________________________________________________
p = polyfit(ratios, a, 2);
if p(1) > 0
    best_temp_grad_ratio = -p(2) / (2 * p(1));
    minT = polyval(p, best_temp_grad_ratio);
else
    [minT, pos] = min(a);
    best_temp_grad_ratio = ratios(pos);
    warning('Quadratic fit opens downward; using data minimum instead.');
end
disp(['Best attraction parameter found (regression): ', num2str(best_temp_grad_ratio)]);
ratios_fit = linspace(min(ratios), max(ratios), 200);
a_fit = polyval(p, ratios_fit);
plot(ratios, a, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Data Points');
hold on;
plot(ratios_fit, a_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Quadratic Fit (Order 2)');
plot(best_temp_grad_ratio, minT, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'DisplayName', 'Calculated Min');
text(best_temp_grad_ratio, minT, ' \leftarrow Minimum', 'FontWeight', 'bold', 'Color', 'r');
title(['Best attraction parameter: ', num2str(best_temp_grad_ratio, '%.4f')]);
xlabel('Ratios');
ylabel('a');
legend('Location', 'best');
grid on;
saveas(gcf, 'quadratic_fit_minimum.png');
pause(5);
close all;
%end of Gemini time________________________________________________________

disp('Full calculation with big image')
starting_image='200x400.bmp';
mkdir('Figure');
mkdir('Topology');
verbose=1;
automate_cell_direct(best_temp_grad_ratio,kp_k0,filling_ratio,heat_sink_temperature,delta_x,p_vol,variation_rate,starting_image,verbose);