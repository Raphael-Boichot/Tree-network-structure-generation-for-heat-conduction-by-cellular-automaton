%https://doi.org/10.1016/j.enconman.2008.09.003
%Raphaël BOICHOT initial code 2009, published 2023
clear;
clc;

kp=[25];%kp/ko ratio
phi=[0.3];%filling ratio

[a]=automate_cell_direct(kp,phi,'400x800.bmp');

