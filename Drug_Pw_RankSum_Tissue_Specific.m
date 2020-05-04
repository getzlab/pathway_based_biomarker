function[] = Drug_Pw_RankSum_Tissue_Specific()

% This function reads pathway activity levels or gene expression levels
% from an excel file alon side with compounds AUC z-score levels and
% calculate the association between pathways or genes and sensitivity to a
% given compound in a specific tiisue type.

fid_out = fopen('Significant_results.txt','wt');
fid_out = fopen('Insignificant_results.txt','wt');

% Data.xlsx id an excel file containing three sheets: (i) microArray - the 
% pathway activity levels across all 460 cell lines calculated from microArray data, 
% (ii)RNAseq - the pathway activity levels across all 460 cell lines
% calculated from RNAseq data, and (iii) AUC_Z - compounds area under the
% curve (AUC) z-score. 
Pw_act = xlsread('DATA.xlsx','microArray','B3:QS811');
[a,Pw_ID] = xlsread('DATA.xlsx','microArray','A3:A811');
[a,Tissue] = xlsread('DATA.xlsx','microArray','B1:QS1');
AUC_Z = xlsread('DATA.xlsx','AUC_Z','B3:QS483');
[a,Compounds] = xlsread('DATA.xlsx','AUC_Z','A3:A483');

U_Tissue = unique(Tissue)';

for kk=1:length(U_Tissue)
    Ind = 1;
    %For every tissue filter all pathways with low variability (i.e
    %pathways with std. and pathway variablity lower than 0.1 
    Act = Pw_act(:,strcmp(U_Tissue(kk), Tissue));
    Act_ID = Pw_ID;
    AUC_TS = AUC_Z(:,strcmp(U_Tissue(kk), Tissue));
    Pval = zeros(20000000,1);
    Results = cell(20000000,2);
    ii=1;
    while (ii < length(Act_ID))
        if(std(Act(ii,:)) < 0.01 || (max(Act(ii,:)) - min(Act(ii,:))) < 0.1)
            Act(ii,:) = [];
            Act_ID(ii) = [];
        else
            ii = ii+1;
        end
    end
    % For every compound - cell lines with AUC z-score under -1.5 are
    % tagged as sensitives and cell lines with z-score above 0 tagged as
    % Not sensitive.
    for ii=1:size(AUC_TS,1)
        for jj=1:size(AUC_TS,2)
            if (AUC_TS(ii,jj) <= -1.5 )
                Group(jj) = {'Sensitive'};
            else
                if(AUC_TS(ii,jj) >= 0 )
                    Group(jj) = {'Not_Sensitive'};
                else
                    Group(jj) = {'NA'};
                end
            end
        end

        if(~isempty(Act(:, strcmp('Sensitive', Group))) && ~isempty(Act(:, strcmp('Not_Sensitive', Group))))
            Sen_G = Act(:, strcmp('Sensitive', Group));
            NotSen_G = Act(:, strcmp('Not_Sensitive', Group));
            
            % In order to perform the statistical test the minimal group
            % size (for either the pathway or the compound) should be 3
            % cell lines.
            if(size(Sen_G,2) >= 3 && size(NotSen_G,2) >= 3)
                for jj=1:length(Act_ID)
                    A = Sen_G(jj,:);
                    B = NotSen_G(jj,:);
                    A(isnan(A)) = [];
                    B(isnan(B)) = [];
                    if(length(A) >= 3 && length(B) >= 3)
                        if(~isnan(ranksum(Sen_G(jj,:), NotSen_G(jj,:),'method','approximate')))
                            Pval(Ind,1) = ranksum(Sen_G(jj,:), NotSen_G(jj,:),'method','approximate');
                            Results(Ind,1) = Compounds(ii);
                            Results(Ind,2) = Act_ID(jj);
                            Ind = Ind+1;
                        end
                    end
                 end    
            end

        end  
    end
    %Visualization of the quantile-quantile plot of the p-values.
    QQ_PLOT(Pval(1:Ind-1));
    %This function calculates the FDR values based on the p-values and
    %print all results with an FDR < 0.25
    Print_significant_results(fid_out,Results(1:Ind-1,:),Pval(1:Ind-1),U_Tissue(kk));
    clear Group Sen_G NotSen_G AUC_TS Act Act_ID Pval Results
end
end
            
        
        

    
    
    

