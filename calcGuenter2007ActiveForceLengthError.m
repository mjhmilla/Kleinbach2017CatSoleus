function meanSqErr = calcGuenter2007ActiveForceLengthError(x,curveSample,errorScaling)

dWdes   = x(1,1);
nuCEdes = x(2,1);
dWasc   = x(3,1);
nuCEasc = x(4,1);

meanSqErr = 0;

for i=1:1:length(curveSample.x)
    flN = 0;
    lceN = curveSample.x(i,1);
    if(curveSample.x(i,1) < 1)
        flN = exp(-abs( (lceN-1)/(dWasc) ).^nuCEasc);   
    else
        flN = exp(-abs( (lceN-1)/(dWdes) ).^nuCEdes);
    end
    err = flN-curveSample.y(i,1);
    meanSqErr = meanSqErr + (err*err);
end
meanSqErr = (meanSqErr ./ length(curveSample.x)) ./ errorScaling;