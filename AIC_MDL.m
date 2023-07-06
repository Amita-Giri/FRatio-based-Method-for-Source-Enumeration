function [AIC_Prediction,MDL_Prediction,AIC,MDL] = AIC_MDL(eigens, p, n, max_sources)
AIC = (zeros(1,max_sources));
MDL = (zeros(1,max_sources));
for k = 0 : max_sources - 1
    prodkp = prod(abs(eigens(k+1:end)).^(1/(p-k)));
    sumkp = sum(abs(eigens(k+1:end)));
    divide = prodkp / ((1/(p-k)) * sumkp);
    logpart = ((p-k)*n) * log(divide);
    linpart = k*(2*p-k);
    AIC(k+1) = -2 * logpart + 2 * linpart;
    MDL(k+1) = -logpart + 0.5 * linpart * log(n);
end
[~,min_aic] = min(AIC);
[~,min_mdl] = min(MDL);
AIC_Prediction = min_aic - 1;
MDL_Prediction = min_mdl - 1;
end