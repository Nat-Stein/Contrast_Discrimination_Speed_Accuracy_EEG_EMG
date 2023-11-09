function [pWilk,pMauchly] = norm_sphere_test(X)
% X is ANOVA table with samples x conditions


% norm_sphere_test
% Test normality and spericity of all stats...


wilk=1;ALPHA=0.05;
clear pWilk
Xconcat = [];
for n = 1:size(X,2)
    for m = 1:size(X,3)
        
        [H, pValue, SWstatistic] = swtest(X(:,n,m), ALPHA,wilk); % Sample X, default alpha = 0.05
        % Shapiro-Wilk test to determine if the null hypothesis of
        %   composite normality is a reasonable assumption regarding the
        %   population distribution of a random sample X.
        % H = 0 => Do not reject the null hypothesis at significance level ALPHA.
        % Inputs:
        %   X - a vector of deviates from an unknown distribution. The observation
        %     number must exceed 3 and less than 5000.
        %
        % Optional inputs:
        %   ALPHA - The significance level for the test (default = 0.05).
        pWilk(n) = pValue;
        hWilk(n) = H;
        Xconcat = [Xconcat X(:,n,m)'];
    end
end
[H, pValue, SWstatistic] = swtest(Xconcat, ALPHA,wilk);
pWilk(n+1) = pValue;
hWilk(n+1) = H;
disp(['HWilk = ' num2str(hWilk)])
disp(['pWilk = ' num2str(pWilk)])



% Mauchly's test of sphericity
[n,p,L,P] = Mauspher(X,ALPHA);
pMauchly = P;
%     Output:
%          n - sample-size.
%          p - variables.
%          L - Mauchly's statistic used to test any deviation from
%              an expected sphericity.
%          P - probability that null Ho: is true.

end