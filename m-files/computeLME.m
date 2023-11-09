function [lmeOut] = computeLME(datset, predVar, fixedE, randomE);
% Computes and displays results of linear mixed-effects model on predVar
% with iterative procedure to test what random slopes significantly improve
% the model fit and then testing whether taking out any fixed effects
% fixedE significantly worsens the model fit
%   Can only take one random effect randomE in its current state

basicLME = [predVar '~'];
for f = 1:length(fixedE)
    if f>1; pls = '+'; else; pls = ''; end
    basicLME = [basicLME pls fixedE{f}];
end
if length(randomE) > 0
    
    for r = 1:length(randomE)
        basicLME = [basicLME '+(1|' randomE{r} ')'];
    end
end
basicMod = fitlme(datset,basicLME);

keepRS = [];
clear altLme
for a = 1:length(fixedE)
    randSlo = {fixedE{a}};
    am = [predVar '~'];
    for f = 1:length(fixedE)
        if f>1; pls = '+'; else; pls = ''; end
        am = [am pls fixedE{f}];
    end
    for r = 1:length(randomE)
        am = [am '+(1'];
        for rs = 1:length(randSlo)
            am = [am '+' randSlo{rs}];
        end
        am = [am '|' randomE{r} ')'];
    end
    
    altLme{a} = am;
    altMod = fitlme(datset,altLme{a});
    results = compare(basicMod,altMod);
    disp([fixedE{a} ' slope: chi^2(' num2str(double(results(2,7))) ') = ' num2str(double(results(2,6))) ', p = ' num2str(double(results(2,8)))]);
    if double(results(2,8)) < 0.05
        keepRS = [keepRS '+' fixedE{a}];
    end
end

finalLME = [predVar '~'];
for f = 1:length(fixedE)
    if f>1; pls = '+'; else; pls = ''; end
    finalLME = [finalLME pls fixedE{f}];
end
if length(randomE) > 0
    for r = 1:length(randomE)
        finalLME = [finalLME '+(1' keepRS '|' randomE{r} ')'];
    end
end
finalMod = fitlme(datset,finalLME);

disp(['Final model: ' finalLME])

clear altLme
for a = 1:length(fixedE)
    am = [predVar '~']; k = 0;
    for f = 1:length(fixedE)
        if f ~= a
            if k>0; pls = '+'; else; pls = ''; end; k = k+1;
            am = [am pls fixedE{f}];
        end
    end
    if length(randomE) > 0
        for r = 1:length(randomE)
            am = [am '+(1' keepRS '|' randomE{r} ')'];
        end
    end
    altLme{a} = am;
    altMod = fitlme(datset,altLme{a});
    results = compare(altMod, finalMod);
    disp([fixedE{a} ': chi^2(' num2str(double(results(2,7))) ') = ' num2str(double(results(2,6))) ', p = ' num2str(double(results(2,8)))]);
    lmeOut{a} = results;
end

end

