function attenuationI = ionAttenuation(layerNum,freq, distance)
    function attenuationConstant = aConstant(layerNum,freq)
    %Calculates the attenuation constant given the ionospheric layer and the frequency of the HF wave.
    
    freq = freq * 10^6;
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

esp = [30 50 70 200];
layerHeight = [80 100 250 350];
phi = atan(2*layerHeight(layerNum)/distance);
kappa = esp(layerNum)/sin(phi);
attenuationI = 10*log(exp(-2*aConstant(layerNum, freq)*kappa*1000))*log10(exp(1));
end

