function meanSqErr = calcGuenter2007PassiveForceLengthError(x,fmax,lopt,curveSample,errorScaling)


nuPee  = x(1,1);
fpee   = x(2,1);
dWdes  = x(3,1);
ellPee0= x(4,1);

lpee0 = ellPee0*lopt;

kpee = (fpee*fmax)/( ( lopt*(dWdes+1-ellPee0) )^nuPee );

meanSqErr = 0;

for i=1:1:length(curveSample.x)
    fpe = 0;
    lceN = curveSample.x(i,1);
    lce = lceN*lopt;
    if(lce >= lpee0)
      fpe = kpee*((lce-lpee0)^nuPee);
    end
    err = (fpe/fmax)-curveSample.y(i,1);
    meanSqErr = meanSqErr + (err*err);
end
meanSqErr = (meanSqErr ./ length(curveSample.x)) ./ errorScaling;