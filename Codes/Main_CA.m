%https://doi.org/10.1016/j.enconman.2008.09.003
%Raphaël BOICHOT initial VB code 2006, initial Matlab code 2009, published 2023
clear;
clc;
close all;

kp_k0=10;%kp/ko ratio, k0=1 by default
phi=0.3;%filling ratio
obj=0;%gradient/temperature ratio for attraction 0=100% gradient, 1=100%temperature 

mkdir('Figure');
mkdir('Topology');
figure('Position',[100 100 800 800]);
[a]=automate_cell_direct(obj,kp_k0,phi,'200x400.bmp');