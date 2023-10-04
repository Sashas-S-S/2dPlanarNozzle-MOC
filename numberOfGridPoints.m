function [gridPoints] = numberOfGridPoints(numberOfMachLines)
    if numberOfMachLines > 1
        gridPoints = (numberOfMachLines-1) + 2 + numberOfGridPoints(numberOfMachLines-1);
    else
        gridPoints = 3;
    end
end