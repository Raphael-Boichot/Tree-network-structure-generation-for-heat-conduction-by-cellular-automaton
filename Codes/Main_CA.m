%https://doi.org/10.1016/j.enconman.2008.09.003
%Raphaël BOICHOT initial VB code 2006, initial Matlab code 2009, published 2023
clear;
clc;
close all;

kp_k0=10;                       %kp/ko ratio, k0=1 by default
filling_ratio=0.3;              %filling ratio
temp_grad_ratio=0.2;            %gradient/temperature ratio for attraction 0=100% gradient, 1=100%temperature 
heat_sink_temperature=298;      %self explanatory
step_x=0.001;                   %step in x/y
p_vol=1e6;                      %surface power or volume density
variation_rate=0.5;             %initial automaton variation rate, decreasing with epoch
starting_image='100x200.bmp';   %self explanatory

mkdir('Figure');
mkdir('Topology');
figure('Position',[100 100 800 800]);
[a]=automate_cell_direct(temp_grad_ratio,kp_k0,filling_ratio,heat_sink_temperature,step_x,p_vol,variation_rate,starting_image);