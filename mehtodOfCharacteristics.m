clear all;
clear;
dThroat = input('Input throat Diameter of Nozzle in mm = ')/1000;

%Chamber conditionsk
Pch = input('Input Chamber Pressure of Nozzle in Pa = ');
Tch = input('Input Chamber Temperature of Nozzle in K = ');

%Ambient Conditions
Pa = input('Input Ambient Pressure of Nozzle in Pa = ');
Ta = input('Input Ambient Temperature of Nozle in K = ');

%% 

%Properties of Gas
gamma = input('Input Cp/Cv of gas = ');
R = 8.314;                                              %J/mol.K
molarMassOfGas = input('Input Molar Mass of gas in g = ')/1000;
Rgas = R/molarMassOfGas;

%%
%Exit Mach
exitMach = input('Input the exit Mach number needed = ');
numberOfMachLines = input('Input number of Centered Expansion Waves = ');

%%

%Throat Conditions
throatMach = 1;                                         %choked flow
throatT = 2*Tch/(gamma+1);
throatP = Pch*((2/(gamma+1))^(gamma/(gamma-1))); 

%total Quantities - constant throughout length of Nozzle
totalP = throatP*(1+((gamma-1)*throatMach^2/2))^(gamma/(gamma-1));
totalT = throatT*(1+((gamma-1)*throatMach^2/2));

gridPoints = ((numberOfMachLines+1)*(numberOfMachLines+2)/2)-1; %numberOfGridPoints(numberOfMachLines) - 1; %theta for first point different for every line

%%
theta = zeros(1,gridPoints);
nu = zeros(1,gridPoints);
K_plus = zeros(1,gridPoints);
K_minus = zeros(1,gridPoints);
mu = zeros(1,gridPoints);
theta_0 = zeros(1,gridPoints);
Ma = zeros(1,gridPoints);

Kmslopes = zeros(1,numberOfMachLines);
Kpslopes = zeros(1,numberOfMachLines);
contourSlopes = zeros(1,numberOfMachLines);

%%
theta_max = prandtlMeyer(exitMach,gamma)/2;
theta_res = theta_max - floor(theta_max);
delta_theta = (theta_max-theta_res)/(numberOfMachLines-1);

bottomWallGridNumber = zeros(1,numberOfMachLines);
topWallGridNumber = zeros(1,numberOfMachLines);

%%
for i=1:numberOfMachLines

bottomWallGridNumber(1,numberOfMachLines-i+1) = gridPoints - (((i+1)*(i+2)/2)-2);
topWallGridNumber(1,numberOfMachLines-i+1) = gridPoints - (((i)*(i+1)/2)-1);

end

xcoords = zeros(1,1+(2*length(bottomWallGridNumber)));
ycoords = zeros(1,1+(2*length(bottomWallGridNumber)));

%%

for i=1:numberOfMachLines
    k=1;
    for j=bottomWallGridNumber(1,i):topWallGridNumber(1,i)
        theta_0(1,j) = theta_res + ((i+k-2)*delta_theta);
        K_minus(1,j) = 2*(theta_0(1,j));
        if(i==1)
            K_plus(1,j) = 0;
            theta(1,j) = (K_plus(1,j) + K_minus(1,j))/2;
            nu(1,j) = (K_minus(1,j)-K_plus(1,j))/2;
            Ma(1,j) = inversePrandtl(nu(1,j),gamma);
            mu(1,j) = asind(1/Ma(1,j));
            if (j==topWallGridNumber(1,i))
                K_minus(1,j) = K_minus(1,j-1);
                K_plus(1,j) = K_plus(1,j-1);
                theta(1,j) = theta(1,j-1);
                nu(1,j) = nu(1,j-1);
                Ma(1,j) = Ma(1,j-1);
                mu(1,j) = mu(1,j-1);
            end
        end
        if(i~=1 && j==bottomWallGridNumber(1,i))
            theta(1,j) = 0;
            nu(1,j) = K_minus(1,j) - theta(1,j);
            K_plus(1,j) = theta(1,j) - nu(1,j);
        end
        if (j~=bottomWallGridNumber(1,i))
            K_plus(1,j) = K_plus(1,j-1);
            theta(1,j) = (K_plus(1,j)+K_minus(1,j))/2;
            nu(1,j) = (K_minus(1,j)-K_plus(1,j))/2;
        end
        Ma(1,j) = inversePrandtl(nu(1,j),gamma);
        mu(1,j) = asind(1/Ma(1,j));
        if(j==topWallGridNumber(1,i) && i~=1)
            K_minus(1,j) = K_minus(1,j-1);
            K_plus(1,j) = K_plus(1,j-1);
            theta(1,j) = theta(1,j-1);
            nu(1,j) = nu(1,j-1);
            Ma(1,j) = Ma(1,j-1);
            mu(1,j) = mu(1,j-1);
        end
        k = k+1;
    end
end

%%
[xcoords,ycoords] = postprocessing(dThroat, theta_max, numberOfMachLines, bottomWallGridNumber, topWallGridNumber, theta, mu, Kpslopes, Kmslopes, contourSlopes, theta_0);

%% Results

exitP = totalP/(1+((gamma-1)*throatMach^2/2))^(gamma/(gamma-1));
exitT = totalT/(1+((gamma-1)*throatMach^2/2));

X = sprintf('Exit Pressure = %d Pa',exitP);
disp(X);
X = sprintf('Exit Temperature = %d K',exitT);
disp(X);

ratio = exitP/Pa;
if ratio<1 && ratio>0.95
    disp('Nozzle is slightly expanded')
else, if ratio<0.95
        disp('Nozzle is over expanded')
else, if ratio > 1
        disp('Nozzle is under expanded')
      end
      end
end   

%% Thrust
chokedMassFlowRate = (Pch*pi*dThroat*dThroat/4)*sqrt(gamma/R/Tch*molarMassOfGas)*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)));
X = sprintf('Mass flow = %d kg/s',chokedMassFlowRate);
disp(X)

Thrust = chokedMassFlowRate*(exitMach) + ((exitP-Pa)*pi*ycoords(1+(2*length(bottomWallGridNumber)))^2/4);
X = sprintf('Thrust = %d N',Thrust);
disp(X)
