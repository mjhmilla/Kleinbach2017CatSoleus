

clc;
close all;
clear all;
fitCrossBridgeStiffnessDampingToKirch199490Hz=0;
flag_useFixedLambdaECM    = 0;

flag_plotEveryCurve       = 0;
flag_plotForceLengthDetail= 0;
flag_makeAndSavePubPlots  = 1;


pubOutputFolder = 'output/plots/MuscleCurves/';
postprocessingDirectoryTree      = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );
addpath('colornames/');


[palette,rgb] = colornames('Alphabet');
colorGuentherSchmittWank2007 = rgb(12,:);

flag_useOctave = 0;
flag_enableNumericallyNonZeroGradients    = 1;

% 1. Active force length curve vs. data
% Solution: There were some initial descrepencies between the experimental force
%length data and a theoretical curve. These errors almost completely go
%away if it is assumed that the experimental recordings are of total 
%path length, rather than fiber length. In this case, when the elastiticy
%of the tendon is taken into account the theoretical active-force-length 
%curve and the transformed data nicely align.

%Failed attempt:
%This creates a cat soleus with an optimal fiber length of 58 mm: this
%is simply way too big to be realistic (given the data I'm working from)
flag_solveForOptimalFiberLengthOfBestFit  = 0; 

%Failed attempt:
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;%...
%  (1/3)*( (1.154-1.087) + (1.23-1.162) + (1.077-1.039) );

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));

%%
% Add the directories needed to run this script
%%
parametersDirectoryTreeMTParams     = genpath('parameters');
parametersDirectoryTreeExperiments  = genpath('experiments');
parametersDirectoryTreeModels       = genpath('models');
parametersDirectoryTreeCurves       = genpath('curves');

addpath(parametersDirectoryTreeMTParams);
addpath(parametersDirectoryTreeExperiments);
addpath(parametersDirectoryTreeModels);
addpath(parametersDirectoryTreeCurves);



scaleOptimalFiberLength      = 1.0; 

scaleMaximumIsometricTension = 1;

if(exist('fitCrossBridgeStiffnessDampingToKirch199490Hz','var')==0)
  fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
end

[felineSoleusMusculotendonProperties, ...
 felineSoleusSarcomereProperties,...
 felineSoleusActiveForceLengthData,...
 felineSoleusPassiveForceLengthData] = createFelineSoleus(...                                          
                                          scaleOptimalFiberLength,...
                                          scaleMaximumIsometricTension,...
                                          fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                                          flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createFelineSoleus(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        flag_useOctave); 
                                        
[felineSoleusNormMuscleCurves,...
 felineSoleusMusculotendonPropertiesUpd,...
 felineSoleusSarcomerePropertiesUpd,...
 activeForceLengthCurveAnnotationPoints,...
 felineSoleusActiveForceLengthDataUpd,...
 felineSoleusPassiveForceLengthDataUpd,...
 forceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      felineSoleusMusculotendonProperties,...
      felineSoleusSarcomereProperties,...
      felineSoleusActiveForceLengthData,...
      felineSoleusPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_useFixedLambdaECM,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);

if(flag_plotEveryCurve==1)
    figH = plotStructOfBezierSplines( felineSoleusNormMuscleCurves,...
                                      'Inverse');                          
    
    %%
    % Note the average offset between the active-force-length curve and
    % the transformed data
    %%
    
    xExp = felineSoleusActiveForceLengthDataUpd(2:end,1);
    yExp = felineSoleusActiveForceLengthDataUpd(2:end,2);
    xCurve = zeros(size(xExp));
    
    for i=1:1:length(xExp)
    xCurve(i,1) = calcBezierFcnXGivenY(yExp(i,1), ...
      felineSoleusNormMuscleCurves.activeForceLengthCurve,... 
      xExp(i,1));
    end                                    
    dx = mean(xCurve-xExp);
    %felineSoleusActiveForceLengthDataUpd(:,1)=...
    %  felineSoleusActiveForceLengthDataUpd(:,1)+dx;
    
    disp('Normalized length offset');
    fprintf('%1.6f lce/lopt \n',felineSoleusActiveForceLengthDataUpd(1,1));
    disp('Average error on the descending limb');
    fprintf('%1.6f lce/lopt \n',dx);
    fprintf('%1.6f mm \n',dx*(felineSoleusMusculotendonPropertiesUpd.optimalFiberLength*1000));
    
    lceNStart = felineSoleusActiveForceLengthDataUpd(1,1);
    save('output/structs/normalizedFiberLengthStartHerzogLeonard2002.mat',...
         'lceNStart');
    
    % dl = felineSoleusActiveForceLengthDataUpd(2:end,1)-1;
    % A  = [dl ones(size(dl))];
    % b  = felineSoleusActiveForceLengthDataUpd(2:end,2);
    % 
    % x     = (A'*A)\(A'*b);
    % y0    = x(2,1);
    % dydx0 = x(1,1);
    % 
    % felineSoleusActiveForceLengthLineBestFit = zeros(length(dl),1);
    % 
    % felineSoleusActiveForceLengthLineBestFit(:,1) = ...
    %   felineSoleusActiveForceLengthDataUpd(2:end,1);
    % dl = felineSoleusActiveForceLengthLineBestFit(:,1)-1;
    % 
    % felineSoleusActiveForceLengthLineBestFit(:,2) = dydx0.*dl + y0;
    % 
    % disp('Active force length line of best fit y=(dydx)*(x-1) + y0');
    % fprintf('dydx: %1.3f\n',dydx0);
    % fprintf('  y0: %1.3f\n',y0);
    
    
    
    %%
    % Plot the derivative of the tendon force length curve on top of the
    % stiffness curve
    %% 
    figure(figH.tendonStiffnessCurve);
    
        curveSample = calcBezierYFcnXCurveSampleVector(...
                        felineSoleusNormMuscleCurves.('tendonForceLengthCurve'), 200,[]);
    
        xmin = min(curveSample.x);
        xmax = max(curveSample.x);
        ymin = min(curveSample.y);
        ymax = max(curveSample.y);
    
    subplot(2,2,1);
      plot(curveSample.x, curveSample.dydx,...
        '--','Color',[1,1,1].*0.5,'LineWidth',2);
      hold on;
    
    %%
    %Plot experimental data over top of the curves where it is available.
    %%
    figure(figH.activeForceLengthCurve);
      subplot(2,2,1);  
      plot(  felineSoleusActiveForceLengthDataUpd(:,1),...
           felineSoleusActiveForceLengthDataUpd(:,2),'xb');
      hold on;          
    
    figure(figH.fiberForceLengthCurve);
      subplot(2,2,1);
      plot(   felineSoleusPassiveForceLengthDataUpd(:,1),...
              felineSoleusPassiveForceLengthDataUpd(:,2),'xb');
      hold on;          
  
end

%%
% Plot the passive force length curve of the 2 segment titin model and
% compare it to the passive force length curve
%%
if(flag_plotForceLengthDetail==1)
    lceN0 = calcBezierFcnXGivenY(0, ...
              felineSoleusNormMuscleCurves.fiberForceLengthCurve);
    lceN1 = calcBezierFcnXGivenY(1, ...
              felineSoleusNormMuscleCurves.fiberForceLengthCurve);
            
    npts = 100;
    lceNSeries = [lceN0:((lceN1-lceN0)/(npts-1)):lceN1]';
    
    fNSeries = zeros(npts,4); % fpe, fecm, fIgp, fPevkIgd,
    lNSeries = zeros(npts,4); % lpe, lecm, lIgp, lPevkIgd,
    
    normLengthIgdFixed = ...
      felineSoleusSarcomerePropertiesUpd.IGDFixedNormLengthAtOptimalFiberLength;
    
    normLengthT12ToZ = ...
      felineSoleusSarcomerePropertiesUpd.ZLineToT12NormLengthAtOptimalFiberLength;
    
    
    for i=1:1:npts
      lceN  = lceNSeries(i,1);
      fpeN  = calcBezierYFcnXDerivative(lceN,...
                felineSoleusNormMuscleCurves.fiberForceLengthCurve,0);
      fecmN = calcBezierYFcnXDerivative(lceN*0.5,...
                felineSoleusNormMuscleCurves.forceLengthECMHalfCurve,0);
              
      lIgpPevkN = lceN*0.5 - normLengthIgdFixed - normLengthT12ToZ; 
              
      [lIgpN, lPevkIgdN, fTiN] = calcSeriesSpringStretch(lIgpPevkN,...
                felineSoleusNormMuscleCurves.forceLengthIgpCurve,...
                felineSoleusNormMuscleCurves.forceLengthIgpInverseCurve, ...
                felineSoleusNormMuscleCurves.forceLengthPevkIgdCurve,...
                felineSoleusNormMuscleCurves.forceLengthPevkIgdInverseCurve);          
    
      fNSeries(i,:) = [fpeN,fecmN,fTiN, fTiN];
      lNSeries(i,:) = [lceN, lceN,lIgpN*2, lPevkIgdN*2];
              
    end
    
    
    fig_forceLength = figure;
      plot(lNSeries(:,1),fNSeries(:,1),'k');
      hold on;
      plot(lNSeries(:,1), fNSeries(:,2)+fNSeries(:,3),'b');
      hold on;
      plot(lNSeries(:,1), fNSeries(:,2),'r');
      hold on;
      plot(felineSoleusPassiveForceLengthDataUpd(:,1),...
           felineSoleusPassiveForceLengthDataUpd(:,2),'xb');
      legend('fpe','fecm+fti','fecm','data');
      xlabel('Norm. Length')
      ylabel('Norm. Force');
      

end

defaultFelineSoleus = struct('musculotendon',...
                            felineSoleusMusculotendonPropertiesUpd,...
                            'sarcomere',...
                            felineSoleusSarcomerePropertiesUpd,...
                            'falData',...
                            felineSoleusActiveForceLengthDataUpd,...
                            'fpeData',...
                            felineSoleusPassiveForceLengthDataUpd,...
                            'curves',...
                            felineSoleusNormMuscleCurves);
                      
save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');                      


%%
% Generate publication quality plots
%%
if(flag_makeAndSavePubPlots==1)
  [fig_pubCurves,subPlotPanel] = ...
      plotMuscleCurves( felineSoleusNormMuscleCurves,...
                      activeForceLengthCurveAnnotationPoints,...
                      felineSoleusMusculotendonPropertiesUpd,...
                      felineSoleusSarcomerePropertiesUpd,...
                      felineSoleusActiveForceLengthDataUpd,...
                      felineSoleusPassiveForceLengthDataUpd,...
                      pubOutputFolder);

   subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
    
   %Zajac tendon parameters
   eIso  = felineSoleusMusculotendonPropertiesUpd.tendonStrainAtOneNormForce;
   fToe  = 2/3;
   fIso  = felineSoleusMusculotendonPropertiesUpd.fiso;
   lsee0 = felineSoleusMusculotendonPropertiesUpd.tendonSlackLength;

   %%
   %Guenther, Schmitt, Wank 2007 tendon   
   %%
   kIsoN=(1.375/eIso);
   kIso=kIsoN*(fIso/lsee0);

   %This tendon model has a constraint that the toe region is exactly
   %half of the change in force between 0 and the linear reference point
   dFsee0  = fToe*fIso;
   dUsee   = (eIso + (fToe*2-1)/kIsoN);
   dUseeL  = (dFsee0/fIso)/kIsoN;
   dUseeNL = (dUsee-dUseeL);
   vsee    = dUseeNL/dUseeL;
   kseeNL  = dFsee0/((dUseeNL*lsee0)^vsee);
   kseeL   = dFsee0/(dUseeL*lsee0);  

   disp('Tendon Parameters (Guenter, Schmitt, Wank (2007))');
   fprintf('%1.7f\t%s\n',lsee0,'lSEE0');
   fprintf('%1.7f\t%s\n',dUseeNL,'dUSEEnll');
   fprintf('%1.7f\t%s\n',dUseeL,'duSEEl');
   fprintf('%1.7f\t%s\n',dFsee0,'dFSEE0');

   npts=100;
   lt      = lsee0 + ([0:(1/(npts-1)):1]').*(dUsee*lsee0);
   ft     = zeros(size(lt));

   idxNL       = find(lt(:,1) < (lsee0+dUseeNL));
   ft(idxNL,1) = kseeNL.*((lt(idxNL,:)-lsee0).^vsee);
   idxL        = find(lt(:,1) >= (lsee0+dUseeNL));
   ft(idxL,1)  = dFsee0 + kseeL.*(lt(idxL,:)-(lsee0+dUseeNL));

   
   plot(lt./(lsee0), ft./(fIso),'Color',colorGuentherSchmittWank2007);
   hold on;

   %%
   %Guenther, Schmitt, Wank 2007 active force length curve 
   %%
   subplot('Position',reshape(subPlotPanel(1,2,:),1,4));

   lopt    = felineSoleusMusculotendonPropertiesUpd.optimalFiberLength;
   dWdes   = 0.14;
   nuCEdes = 3.0;
   dWasc   = 0.57;
   nuCEasc = 4.0;
   
   curveSampleFL = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('activeForceLengthCurve'), ...
                  50,felineSoleusNormMuscleCurves.('activeForceLengthCurve').xEnd);
    
   x0           = [dWdes;nuCEdes;dWasc;nuCEasc];
   errorScalingFL = calcGuenter2007ActiveForceLengthError(x0,curveSampleFL,1);
   errFcnFL = @(argX)calcGuenter2007ActiveForceLengthError(argX,curveSampleFL,errorScalingFL);

   lb = zeros(size(x0));
   ub = ones(size(x0)).*100;

   [x,fval,exitflag,output]=fmincon(errFcnFL,x0,[],[],[],[],lb,ub);

   dWdes   = x(1,1);
   nuCEdes = x(2,1);
   dWasc   = x(3,1);
   nuCEasc = x(4,1);   

   disp('Active Force Length Parameters (Guenter, Schmitt, Wank (2007))');
   fprintf('%1.7f\t%s\n',dWdes,'dWdes');
   fprintf('%1.7f\t%s\n',nuCEdes,'nuCEdes');
   fprintf('%1.7f\t%s\n',dWasc,'dWasc');
   fprintf('%1.7f\t%s\n',nuCEasc,'nuCEasc');

   npts=100;
   lce = ([0.5:((1.6-0.5)/(npts-1)):1.6]').*lopt;
   flN  = zeros(size(lce));

   idxAsc = find(lce <= lopt);
   idxDes = find(lce > lopt);
   flN(idxAsc,1) = exp(-abs( ((lce(idxAsc,1)./lopt)-1)/(dWasc) ).^nuCEasc);
   flN(idxDes,1) = exp(-abs( ((lce(idxDes,1)./lopt)-1)/(dWdes) ).^nuCEdes);

   plot(lce./lopt, flN,'Color',colorGuentherSchmittWank2007);
   hold on;

   %%
   %Guenther, Schmitt, Wank 2007 passive force length curve 
   %%
   subplot('Position',reshape(subPlotPanel(1,2,:),1,4));

   lopt    = felineSoleusMusculotendonPropertiesUpd.optimalFiberLength;   
   fIso  = felineSoleusMusculotendonPropertiesUpd.fiso;

   ellPee0 = forceLengthCurveSettings.normLengthZero;
   nuPee = 2.5;
   fpee  = 1.0;

   dWdes = 0.6;

   x0 = [nuPee;fpee;dWdes;ellPee0];
   lb = zeros(size(x0));
   ub = ones(size(x0)).*50;

   curveSampleFPE = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('fiberForceLengthCurve'), ...
                  50,felineSoleusNormMuscleCurves.('fiberForceLengthCurve').xEnd);
    
   errorScalingFPE = calcGuenter2007PassiveForceLengthError(x0,fIso,...
       lopt,curveSampleFPE,1);
   errFcnFPE = @(argX)calcGuenter2007PassiveForceLengthError(argX,fIso,...
       lopt,curveSampleFPE,errorScalingFPE);

   [x,fval,exitflag,output]=fmincon(errFcnFPE,x0,[],[],[],[],lb,ub);

   nuPee  = x(1,1);
   fpee   = x(2,1);
   dWdes  = x(3,1);
   ellPee0= x(4,1);

   lpee0 = ellPee0*lopt;   
   kpee = (fpee*fIso)/( ( lopt*(dWdes+1-ellPee0) )^nuPee );

   npts=100;
   lceN = [ellPee0*0.5:(1.4-ellPee0*0.5)/(npts-1):1.4]';
   lce  = lceN.*lopt;
   fpe  = zeros(size(lce));

   for i=1:1:length(lce)
        if(lce(i,1) >= lpee0)
            fpe(i,1) = kpee*((lce(i,1)-lpee0)^nuPee);
        end
   end

   plot(lceN,fpe./fIso,'Color',colorGuentherSchmittWank2007);
   hold on;

   disp('Parallel Elastic Element Parameters (Guenter, Schmitt, Wank (2007))');
   fprintf('%1.7f\t%s\n',nuPee,'nuPee');
   fprintf('%1.7f\t%s\n',fpee,'fpee');
   fprintf('%1.7f\t%s\n',dWdes,'dWdes');
   fprintf('%1.7f\t%s\n',ellPee0,'ellPee0');

   %%
   %Guenther, Schmitt, Wank 2007 force-velocity curve 
   %%
   subplot('Position',reshape(subPlotPanel(1,3,:),1,4));

   vmaxN = felineSoleusMusculotendonPropertiesUpd.maximumNormalizedFiberVelocity;

    Arel0=0.1;
    Brel0=1.0;
    Seec=2.0;
    Feec=1.8;

    x0 = [Arel0;Brel0;Seec;Feec];
    lb = zeros(size(x0));
    ub = ones(size(x0)).*50;

    curveSampleFv = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('fiberForceVelocityCurve'), ...
                  50,felineSoleusNormMuscleCurves.('fiberForceVelocityCurve').xEnd);

    errorScalingFv = calcGuenter2007ForceVelocityError(x0,fIso,...
       lopt,vmaxN,curveSampleFv,1);
    errFcnFv = @(argX)calcGuenter2007ForceVelocityError(argX,fIso,...
       lopt,vmaxN,curveSampleFv,errorScalingFv);

   [x,fval,exitflag,output]=fmincon(errFcnFv,x0,[],[],[],[],lb,ub);

   Arel0 = x(1,1);
   Brel0 = x(2,1);
   Seec  = x(3,1);
   Feec  = x(4,1);

   disp('Force Velocity Parameters (Guenter, Schmitt, Wank (2007))');
   fprintf('%1.7f\t%s\n',Arel0,'Arel0');
   fprintf('%1.7f\t%s\n',Brel0,'Brel0');
   fprintf('%1.7f\t%s\n',Seec,'Seec');
   fprintf('%1.7f\t%s\n',Feec,'Feec');
   fprintf('%1.7f\t%s\n',vmaxN,'vmaxN');


   %Constants
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
   
   %Sample the curve
   vceN = curveSampleFv.x;
   vce  = vceN.*(lopt*vmaxN);
   fv   = zeros(size(vce));

   for i=1:1:length(fv) 
    if(vce(i,1) <= 0)
      fv(i,1) = fIso*( (q*flN + Arel) / (1 - (vce(i,1)/(Brel*lopt)) )  - Arel );
    else
      fv(i,1) = fIso*( (q*flN + Arele) / (1 - (vce(i,1)/(Brele*lopt)) )  - Arele );
    end
   end

   plot(vceN,fv./fIso,'Color',colorGuentherSchmittWank2007);
   hold on;

   print('-dpdf', [pubOutputFolder,'fig_Pub_MuscleCurves.pdf']); 
   here=1;
end

   
%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
