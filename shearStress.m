function [C_shear] = shearStress(alpha, rdiv)

global rim mat tauArr

a = zeros(1,length(rim));
g = zeros(1,length(rim));
C = zeros(1,length(rim));

for k = 1:length(rim)-1
    g(k) = mat.Q{k}(4,4);
    a(k) = (mat.rho{k} * alpha) / g(k);
end

for k = length(rim):-1:2
    C(k-1) = (a(k-1) - (g(k)*a(k)/g(k-1))) * (rim(k)^4/8) + (g(k)*C(k))/g(k-1);
end
C_shear = C;

for k = 1:length(rim) - 1 
    dr = linspace(rim(k),rim(k+1),rdiv);
    rStart = (k-1)*rdiv + 1;
    rEnd = k*rdiv;
    tauArr(1,rStart:rEnd) = (-mat.Q{1,k}(4,4) * ((a(k)*dr.^2)./4 - (2*C(k))./dr.^2));
end

