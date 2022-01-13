function meanSqErr = calcGuenter2007ForceVelocityError(x,fmax,lopt,vmaxN,curveSample,errorScaling)


Arel0 = x(1,1);
Brel0 = x(2,1);
Seec  = x(3,1);
Feec  = x(4,1);

q     = 1;
lceN  = 1; 
flN   = 1;

LArel = 1;
if(lceN > 1)
    LArel = flN;
end

LBrel = 1;
QArel = (1/4)*(1+3*q);
QBrel = (1/7)*(3+4*q);

Arel = Arel0*LArel*QArel;
Brel = Brel0*LBrel*QBrel;

Arele = -Feec*q*flN;
Brele = Brel*(1-Feec) / (Seec*(1 + (Arel)/(q*flN)));

meanSqErr = 0;

for i=1:1:length(curveSample.x)
    fv = 0;
    vceN = curveSample.x(i,1);
    vce = vceN*lopt*vmaxN;
    if(vce <= 0)
      fv = fmax*( (q*flN + Arel) / (1 - (vce/(Brel*lopt)) )  - Arel );
    else
      fv = fmax*( (q*flN + Arele) / (1 - (vce/(Brele*lopt)) )  - Arele );
    end
    fvN = fv/fmax;
    err = (fvN)-curveSample.y(i,1);
    meanSqErr = meanSqErr + (err*err);
end
meanSqErr = (meanSqErr ./ length(curveSample.x)) ./ errorScaling;