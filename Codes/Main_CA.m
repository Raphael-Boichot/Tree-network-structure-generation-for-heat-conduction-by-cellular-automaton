%https://doi.org/10.1016/j.enconman.2008.09.003
%Raphaël BOICHOT initial VB code 2006, initial Matlab code 2009, published 2023
clear;
clc;

kp=[100];%kp/ko ratio
phi=[0.3];%filling ratio

[a]=automate_cell_direct(kp,phi,'100x200.bmp');

