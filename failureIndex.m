function [strRatio, peakStr, peakLoc] = failureIndex(rdiv,b)

global results rArr rim mat


for k = 1:length(rim) - 1
    rStart = (k-1)*rdiv + 1;
    rEnd = k*rdiv;

    LgTens = mat.stren{k}(1); % longitudinal tensile strength
    LgComp = mat.stren{k}(2); % longitudinal compressive strength
    TranTens = mat.stren{k}(3); % transverse tensile strength
    TranComp = mat.stren{k}(4); % transverse compressive strength
    sh6 = mat.stren{k}(5); % longitudinal shear strength
%     sh4 = mat.stren{k}(6); % transverse shear

    F11 = 1/(LgTens*LgComp); % F11
    F1 = 1/LgTens - 1/LgComp; % F1
    F33 = 1/(TranTens*TranComp); % F33
    F3 = 1/TranTens - 1/TranComp; % F3
    F22 = F33;
    F2 = F3;
    F13 = -0.5*sqrt(F11*F33); % F13
    F12 = F13;
    F66 = 1/sh6^2; % F66
    F23 = F66; %F22 - 1/(2*sh6^2);
    

    sigt = results.sArr{b}(1,rStart:rEnd); % sig1 circumferential stress
    sigax = results.sArr{b}(2,rStart:rEnd); % sig2 axial stress
    sigr = results.sArr{b}(3,rStart:rEnd); % sig3 radial
    %     tau = results.tauArr{b}(1,rStart:rEnd); % tau12

    A = F11*sigt.^2 + F22*sigax.^2 + F33*sigr.^2 +...
        2*F13*sigt.*sigr + 2*F12.*sigt.*sigax + 2*F23.*sigax.*sigr;% + F66*tau.^2;
    B = F1*sigt + F2*sigax + F3*sigr;
    C = -1;
    
    R(1,rStart:rEnd) = (-B + sqrt(B.^2 - 4.*A.*C)) ./ (2*A);
    strRatio(1,rStart:rEnd) = R(1,rStart:rEnd).^-1;

end
[peakStr(b), ind] = max(strRatio(1,:));
peakLoc(b) = rArr(ind);
