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