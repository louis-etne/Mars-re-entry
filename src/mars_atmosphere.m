clear all;
close all;
clc;

% NASA model stops at 100km
h1 = 0:2:100;
% Teacher models stop at 150km but is less precise
h2 = 110:10:150;
% So we use the NASA model first and after the teacher
h = [h1 h2] * 1000; % Convert h in meters

% Temperature is in K
t1 = [214, 213.8, 213.4, 212.4, 209.2, 205, 201.4, 197.8,...
      194.6, 191.4, 188.2, 185.2, 182.5, 180, 177.5, 175, 172.5, 170,...
      167.5, 164.8, 162.4, 160, 158, 156, 154.1, 152.2, 150.3, 148.7,...
      147.2, 145.7, 144.2, 143, 142, 141, 140, 139.5, 139, 139, 139,...
      139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139];
t2 = [149.4, 159.7, 170, 245.1, 288.6];
T = [t1 t2];