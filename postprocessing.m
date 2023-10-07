function [xcoords,ycoords] = postprocessing(dThroat, theta_max, numberOfMachLines, bottomWallGridNumber, topWallGridNumber, theta, mu, Kpslopes, Kmslopes, contourSlopes, theta_0)

    for i=1:length(bottomWallGridNumber)
        Kpslopes(1,i) = (theta(1,bottomWallGridNumber(1,i)) + theta(1,topWallGridNumber(1,i)) + mu(1,bottomWallGridNumber(1,i)) + mu(1,topWallGridNumber(1,i)))/2;
        Kmslopes(1,i) = theta_0(1,i) - 90;
        if(i==1)
            contourSlopes(1,i) = (theta_max + theta(1,topWallGridNumber(1,i)))/2;
        else
            contourSlopes(1,i) = (theta(1,topWallGridNumber(1,i-1))+theta(1,topWallGridNumber(1,i)))/2;
        end
    end

    xcoords(1,1) = 0;
    ycoords(1,1) = dThroat/2;
    
    %xcoords(1,1) = 0;
    %ycoords(1,1) = dThroat/2;
    
    %%
    for i=1:numberOfMachLines
        ycoords(1,i+1) = 0;
        xcoords(1,i+1) = -ycoords(1,1)/tand(Kmslopes(1,i));
        if (i==1)
            A = [-tand(contourSlopes(i)) 1; -tand(Kpslopes(i)) 1];
            B = [ycoords(1); -tand(Kpslopes(i))*xcoords(i+1)];
            X = linsolve(A,B);
            xcoords(i+1+length(bottomWallGridNumber)) = X(1);
            ycoords(i+1+length(bottomWallGridNumber)) = X(2);
            continue;
        end
        A = [-tand(contourSlopes(i)) 1; -tand(Kpslopes(i)) 1];
        B = [(-tand(contourSlopes(i))*xcoords(i+length(bottomWallGridNumber)))+ycoords(i+length(bottomWallGridNumber)); -tand(Kpslopes(i))*xcoords(i+1)];
        X = linsolve(A,B);
        xcoords(i+1+length(bottomWallGridNumber)) = X(1);
        ycoords(i+1+length(bottomWallGridNumber)) = X(2);
    end
    x = [0 0];
    y = [0 0];
    %%
    hold off;
    figure(1)
    for i=1:numberOfMachLines
        y(1) = ycoords(1);
        y(2) = ycoords(i+1);
        x(1) = xcoords(1);
        x(2) = xcoords(i+1);
        plot(x,y)
        hold on;
    end

%%
    for i=1:numberOfMachLines
        y(1) = ycoords(i+1);
        y(2) = ycoords(i+1+length(bottomWallGridNumber));
        x(1) = xcoords(i+1);
        x(2) = xcoords(i+1+length(bottomWallGridNumber));
        plot(x,y)
        hold on;
    end
end