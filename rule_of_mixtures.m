ff = .665; %fiber volume fraction
mf = 1-ff; %matrix volume fraction

ef = 310; %fiber tensile modulus [gpa]
gf = 84; %fiber shear modulus [gpa]
sf = 5.3; %fiber tensile strength [gpa]
vf = 0.2;

em = 2.1; %matrix elastic modulus [gpa]
gm = em / sqrt(3); %matrix shear modulus based on Jordan et al. 2008 (Mech. Prop. Epon826/DEA)
sm = .075; %matrix tensile strength [gpa]
vm = 0.38; % poissons ratio

e11 = (ff * ef) + mf * em
e22 = ((ff / ef) + (mf / em))^-1
sh = (gf * gm) / ((gf*vm) + (gm*vf))

x11t = sf * ff + sm * mf
x11c = ((ff / sf) + (mf / sm))^-1

x22t = (1 - (sqrt(ff) - ff) * (1-em/ef)) * sm
tau = gm/mf