function[] = Down_sampling()

[Pw_act, info] = xlsread('E:\PostDoc\Pw based identification of SNPs\Results_Excel\New Analysis new PathOlogist.xlsx','microArray');
Pw_ID = info(3:size(info,1),1);
Tissue = info(1,2:size(info,2))';
Cells = info(2,2:size(info,2))';
[AUC_Z, info] = xlsread('E:\PostDoc\Pw based identification of SNPs\Results_Excel\New Analysis new PathOlogist.xlsx','AUC_Z');
Compounds = info(3:size(info,1),1);
clear info
[Row, Col] = size(AUC_Z); 
U_Tissue = unique(Tissue);

%Choose specific tissue
NSCLC_Act = Pw_act(:,strcmp(U_Tissue(8), Tissue));
NSCLC_AUC = AUC_Z(:,strcmp(U_Tissue(8), Tissue));

% Classify cell lines in every compound
for ii = 1:size(NSCLC_AUC,1)
    for jj = 1:size(NSCLC_AUC,2)
        if (NSCLC_AUC(ii, jj) <= -1.5 )
            Group(ii, jj) = {'Sensitive'};
        else
            if(NSCLC_AUC(ii, jj) >= 0 )
                Group(ii, jj) = {'Not_Sensitive'};
            else
                Group(ii, jj) = {'NA'};
            end
        end
    end
end

C = 1;
while (C < length(Pw_ID))
    if(std(NSCLC_Act(C,:)) <= 0.01 || (max(NSCLC_Act(C,:)) - min(NSCLC_Act(C,:))) <= 0.1 )
        NSCLC_Act(C,:) = [];
        Pw_ID(C) = [];
    else
        C = C+1;
    end
end

Count = 1;
Results = zeros(4120,2);

for ii = 10:size(NSCLC_Act,2) % Sample size
    for jj=1:10 % number of iterations
        Pval = zeros(5000000,1);
        P_ind = 1;
        Rand_p = randperm(size(NSCLC_Act,2), ii);
        ACT = NSCLC_Act(:,Rand_p);
        AUC = Group(:,Rand_p);  
        ID = Pw_ID;
        
        for kk = 1:size(AUC,1)
            if(~isempty(ACT(:, strcmp('Sensitive', AUC(kk,:)))) && ~isempty(ACT(:, strcmp('Not_Sensitive', AUC(kk,:)))))
                Sen_G = ACT(:, strcmp('Sensitive', AUC(kk,:)));
                NotSen_G = ACT(:, strcmp('Not_Sensitive', AUC(kk,:)));
            
                if(size(Sen_G,2) >= 3 && size(NotSen_G,2) >= 3)
                    for tt=1:length(ID)
                        Pval(P_ind,1) = ranksum(Sen_G(tt,:), NotSen_G(tt,:),'method','approximate');
                        P_ind = P_ind+1;
                     end    
                end

            end  
        end
        
        Pval = Pval(1:P_ind-1,:);
        Pval(isnan(Pval)) = [];
        Pval(:,2) = mafdr(Pval(:,1), 'BHFDR', 'TRUE');
        Results(Count,1) = ii;
        Results(Count,2) = length(find(Pval(:,2) < 0.25));
        Count = Count+1;
        clear Sen_G NotSen_G AUC ACT ID Pval Rand_p 
    end       
end

Results = Results(1:Count,:);
end
            
        
        

    
    
    

