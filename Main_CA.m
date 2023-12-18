%https://doi.org/10.1016/j.enconman.2008.09.003
%Raphaël BOICHOT initial VB code 2006, initial Matlab code 2009, published 2023
clear;
clc;

kp=10;%kp/ko ratio
phi=0.3;%filling ratio
obj=0.05;%gradient/temperature ratio for attraction 0=100% gradient, 1=100%temperature 

mkdir('Figure');
mkdir('Topology');

[a]=automate_cell_direct(obj,kp,phi,'200x400.bmp');

