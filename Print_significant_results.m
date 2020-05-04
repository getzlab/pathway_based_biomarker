function[] = Print_significant_results(fid_out,Info,Pval,Source)

FDR = mafdr(Pval, 'BHFDR','TRUE');
Loc = find(FDR<0.25);
for ii=1:length(Loc)
    fprintf(fid_out,'%s\t%s\t%s\t%d\n',cell2mat(Source), cell2mat(Info(Loc(ii),1)), cell2mat(Info(Loc(ii),2)),FDR(Loc(ii)));
end

end