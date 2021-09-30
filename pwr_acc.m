rim = [.16 .2 .33];
w = 24150 * pi / 30;
p = linspace(100, 100e6, 100);
rho = [2810, 1560];
h = .43;

for k = 1:length(rim) - 1
    in(k) = 0.5 * rho(k) * pi * h * (rim(k+1)^4 - rim(k)^4);
end
I = sum(in);

a = p / (I * w);

plot(p,a*(60^2/(2*pi)))