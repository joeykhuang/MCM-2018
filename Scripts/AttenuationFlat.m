function attenuationGF = AttenuationFlat(phi,permittivity,resistivity,freq,fOrR,sOrG)
    terrainRuggednessP = 20;
    lambda = 3*10^8 / (freq*10^6);
        function polarization = Polarization(hOrV,phi,permittivity,conductivity,freq)
        wavelength = 3*10^8 / (freq * 10^6);
        permPrime = permittivity - 60* wavelength * conductivity * 1i;
        srRoot = sqrt(permPrime - (cos(phi))^2);
        %Returns the polarization of the HF wave.
            if hOrV == 'h'
                polarization = (sin(phi) - srRoot)/(sin(phi) + srRoot);
            elseif hOrV == 'v'
                polarization = (permPrime* sin(phi) - srRoot)/(permPrime*sin(phi) + srRoot);
            end
        end
    hPol = abs(Polarization('h',phi,permittivity,resistivity,freq));
    vPol = abs(Polarization('v',phi,permittivity,resistivity,freq));
    if fOrR == 'f'
        attenuationGF = 10*log10(((abs(hPol)^2+abs(vPol)^2))/2);
    elseif sOrG == 's' && fOrR == 'r'
        attenuationGF = - 4;
    elseif sOrG == 'g' && fOrR == 'r'
        atGFH = 10*log10(abs(hPol)^2*exp(-(terrainRuggednessP * sin(phi)/ lambda)^2));
        atGFV = 10*log10(abs(vPol)^2*exp(-(terrainRuggednessP * sin(phi)/ lambda)^2));
        attenuationGF = (atGFH + atGFV )/ 2; 
    end
end