function [s] = IM7_8552_Tzeng2001(mstiff)
% This function reproduces the compliance curve for IM7/8552 discussed by
% Tzeng 2001 in "Viscoelastic Analysis of Composite Cylinders Subjected to
% Rotation." The compliance was as a best fit equation in eq.19 of this
% paper. 
%
% The units given in the paper are inches and pounds. The compliance is
% in^2/lb

global t

s(1) = .00040679;
s(2) = .00519368 * t^0.03;
s(3) = s(2);
s(4) = 0.0095382 * t^0.03;
s(5) = 0.3;
s(6) = 0.36;
