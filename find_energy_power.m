function [E, P] = find_energy_power(h, tStep, alpha, accType)
global rim w mat

in = zeros(1, length(rim)-1);
en = zeros(2, length(rim) -1);

if strcmp(accType, 'const')
     w2 = w + alpha(tStep);
 else
     %under construcrtion
end
 
for b = 1:length(rim) - 1
    in(b) = 0.5 * mat.rho{b} * pi * h * (rim(b+1)^4 - rim(b)^4);
    en(1,b) = 0.5 * in(b) * w^2;
    en(2,b) = 0.5 * in(b) * w2^2;
end
I = sum(in);
E = sum(en, 2);
P = (E(2) - E(1)) / tStep;
 
 
     
 