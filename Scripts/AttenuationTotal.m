function S_R = AttenuationTotal(numOfHops,layer,distance,permittivity,resistivity,initPower,fOrR,sOrG,oneHop, SSN, month, time, transLat, transLon, recLat, recLon, dayTime, sunTime)
%Calculates the SNR and the total attenuation of a HF wave
    function MUFVal = MUF(SSN, month, time, transLat, transLon, recLat, recLon, dayTime, sunTime)
    %MUF calculation
        if 6 < time && time < 18
            Isday = 'd';
        else
            Isday = 'n';
        end
        function A2Return = A2(SSN)
            A2Return = 1.3022 - (0.00156)*SSN;
        end
        function A3Return = A3(month)
            m = month*2*pi/12;
            A3Return = 0.9925+0.011*sin(m) + 0.087*cos(m) + 0.043*sin(2*m)+ 0.003*cos(2*m) - 0.013*sin(3*m)...
                               -0.022*cos(3*m) + 0.003*sin(4*m) + 0.005*sin(5*m) + 0.018*cos(6*m);
        end
        function A4Return = A4(dOrN, localtime,sunsetTime)
            if dOrN == 'd'
                A4Return = 1.11 + 0.01 * localtime;
            elseif dOrN == 'n'
                t = localtime - sunsetTime;
                A4Return = 1.0195 - 0.06*sin(2*t) - 0.037*cos(2*t) + 0.018*sin(4*t) - 0.003*cos(4*t) ...
                            + 0.025*sin(6*t) + 0.018*cos(6*t) + 0.007*sin(8*t) + 0.005*cos(8*t) + 0.006*sin(10*t)...
                            + 0.017*cos(10*t) - 0.009*sin(12*t) - 0.004*cos(12*t);
            end
        end
        function MReturn = M(transLat, transLon, recLat, recLon)
            function greatCirDis = Haversine(transLat, transLon, recLat, recLon)
                phiT = degtorad(transLat);
                phiR = degtorad(recLat);
                lambdaT = degtorad(transLon);
                lambdaR = degtorad(recLon);
                greatCirDis = 2*asin(sqrt((sin((phiR - phiT)/2))^2+cos(phiT)*cos(phiR)*(sin(((lambdaR-lambdaT)/2)))^2));
            end
            d = Haversine(transLat, transLon, recLat, recLon)*6371;
            if d <= 4000
                k = 1;
            elseif d > 4000
                k = 0.5;
            end
            MReturn = (1+2.5*(sin(2.5*Haversine(transLat, transLon, recLat, recLon)*k))^1.5)*(1 - 0.1*exp((dayTime - 24)/3));
        end
        function f0F2Return = f0F2(SSN)
            f0F2Return = sqrt(2+(0.814*SSN+22.23)*sqrt(0.5));
        end

        MUFVal = A2(SSN)*A3(month)*A4(Isday,time,sunTime)*M(transLat, transLon, recLat, recLon)*f0F2(SSN);
    end
    
    %Attenuations functions
    function AttenuationFS = FSPL(distance,freq)
    %Calculaties the Free Space Path Loss for a HF wave.
    AttenuationFS = 20*log10(freq) + 20*log10(numOfHops*distance) + 32.45;
    end

    function attenuationI = ionAttenuation(layerNum)
        function attenuationConstant = aConstant(layerNum)
        %Calculates the ionospheric attenuation given the reflection layer of the HF wave.

        colliFreq = [1000000,10000,1000,1000];
        eDensity = [10^10, 10^11, 10^12, 10^12];
        theta = colliFreq(layerNum);
        electronDensity = eDensity(layerNum);
        omegap = sqrt((electronCharge ^ 2 * electronDensity)/(electronMass * epsilonZero));
        attenuationConstant = (omegap.^2*theta)/(2 * 3 * 10^8 * (omega.^2 - theta ^ 2));
        end

    esp = [30 50 70 200];
    kappa = esp(layerNum)/sin(phi);
    attenuationI = 10*log(exp(-2*aConstant(layerNum)*kappa*1000))*log10(exp(1));
    end

    function attenuationGF = AttenuationFlat(phi,permittivity,resistivity,freq,fOrR,sOrG)
    %Calculates the ground attenuation
    terrainRuggednessP = 20;
    lambda = 3*10^8 / (freq*10^6);
        function polarization = Polarization(hOrV,phi,permittivity,resistivity,freq)
        wavelength = 3*10^8 / (freq * 10^6);
        permPrime = permittivity - 60* wavelength * resistivity * 1i;
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
    
    function S_R = SNR(iniPower,totalA,noiseLevel)
        iniPower = 10*log10(iniPower);
        S_R = iniPower + totalA + 30 + noiseLevel;
    end
    %Attenuations functions END

    function skipZone = SkipZone(MUF, critFreq, tlayerHeight)
        skipZone = 2*tlayerHeight*sqrt((MUF/critFreq)^2 - 1);
    end

    %Constants
    if layer == 'D'
        layerNum = 1;
    elseif layer == 'E'
        layerNum = 2;
    elseif layer == 'F1'
        layerNum = 3;
    elseif layer == 'F2'
        layerNum = 4;
    end
    
    
    MUFVal = MUF(SSN, month, time, transLat, transLon, recLat, recLon, dayTime, sunTime);
    frequency = 0.85*MUFVal;
    layerHeight = [80 100 250 350];
    phi = atan(2*layerHeight(layerNum)/(distance));
    omega = 2 * pi .* frequency * 10^6;
    electronCharge = 1.6 * 10^-(19);
    electronMass = 9 * 10 ^-31;
    epsilonZero = 8.854 * 10 ^ -12;
    terrainRuggednessP = 20;
    %Constants END
    if distance < SkipZone(MUFVal, 5, 30)
        totalAttenuation = 0; 
    else
        %Calculations 
        aIArray = [0 0 0 0];
        for layerNums = 1:4
            aIArray(layerNums) = aIArray(layerNums) + ionAttenuation(layerNums);
        end    

        A_fs = FSPL(distance * numOfHops, frequency);
        A_i = 0;
        for layerProg = 1:layerNum - 1
            A_i = A_i + aIArray(layerProg);
        end
        A_g = AttenuationFlat(phi,permittivity,resistivity,frequency,fOrR,sOrG);
        if fOrR == 'f' || fOrR == 'r' && oneHop == 'n'
            totalAttenuation = -A_fs + numOfHops*A_i + (numOfHops)*A_g - 7;
        elseif fOrR == 'r' && oneHop == 'y'
            totalAttenuation = -A_fs + numOfHops*A_i + A_g + (numOfHops - 1)*AttenuationFlat(phi,permittivity,resistivity,frequency,'f',sOrG) - 7;
        end
        S_R = SNR(initPower,totalAttenuation,139);
        receivedStrength = 10*log10(initPower) + totalAttenuation;
        %Calculations END      
   end      
   
end

