function [varargout] = IM7_8552_Saleeb2003(mstiff)
% This fuction produces the time dependent compliance IM7/8552 composite studied
% by Saleeb et. al in "A study of time-dependent and anisotropic effects on the
% deformation response of two flywheel designs." The compliance was caluclated 
% from the strain data in figure 2f useing a pdf to graph converter and curve 
% fitting the discrete data using matlab's built in curve fitting toolbox. 
% The data was fit uings a two term exponential function the goodness of fit 
% parameters can be found below. The raw data can be found in 
% 'IM7 creep strain.csv'
% GOF: adj-r^2: 0.9895, SSE: 4.5e-24, RMSE: 5.5e-13
% This function uses the total time to predict the compliance. This is
% taken as the S22 = S33, compliance as of writing this. The raw data is
% curve fit to linear time so t should be linear to a maximum of 252518
% seconds.

global t

a = 1.395e-10;
b = 2.912e-7;
c = -7.472e-12;
d = -3.708e-5;

s(1) = mstiff(1);
s(2) = a.*exp(b.*t) + c.*exp(d.*t);
s(3) = s(2);
s(4) = mstiff(4);
s(5) = mstiff(5);
s(6) = s(5);
varargout{1} = s;
