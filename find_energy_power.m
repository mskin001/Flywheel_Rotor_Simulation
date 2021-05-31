function [E] = find_energy_power(h)
global rim w mat

in = zeros(1, length(rim)-1);
en = zeros(2, length(rim) -1);
 
for b = 1:length(rim) - 1
    in(b) = 0.5 * mat.rho{b} * pi * h * (rim(b+1)^4 - rim(b)^4);
    en(b) = 0.5 * in(b) * w^2;
end
I = sum(in);
E = sum(en);
 
 
     
 