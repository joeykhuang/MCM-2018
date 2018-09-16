function attenuationConstant = aConstant(layerNum,freq)
%Calculates the attenuation constant given the ionospheric layer and the frequency of the HF wave.

colliFreq = [1000000,10000,1000,1000];
eDensity = [10^10, 10^11, 10^12, 10^12];
theta = colliFreq(layerNum);
electronDensity = eDensity(layerNum);
omega = 2 * pi * freq;
electronCharge = 1.6 * 10^-(19);
electronMass = 9 * 10 ^-31;
epsilonZero = 8.854 * 10 ^ -12;
omegap = sqrt((electronCharge ^ 2 * electronDensity)/(electronMass * epsilonZero));
attenuationConstant = (omegap^2*theta)/(2 * 3 * 10^8 * (omega^2 - theta ^ 2));
end

