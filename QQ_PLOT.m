function[] = QQ_PLOT(Pval)

Observed = sort(-log10(Pval));
for ii=1:length(Observed)
    Expected(ii) = -log10(ii/length(Observed));
end
Expected = sort(Expected);
h = qqplot(Expected, Observed);       
end
