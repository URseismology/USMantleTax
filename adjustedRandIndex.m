function ARI = adjustedRandIndex(true_labels, predicted_labels)
    % Create a contingency table
    [C,~,~] = contingency(true_labels, predicted_labels);
    
    % Sum over rows & columns
    nis = sum(C,2);
    njs = sum(C,1);
    
    % Sum over pairs of bins
    nij = sum(sum(C.*(C-1)/2));
    
    % Adjusted for chance indices
    ni = sum(nis.*(nis-1)/2);
    nj = sum(njs.*(njs-1)/2);
    
    % The number of all pairs
    n = sum(nis)*((sum(nis)-1)/2);
    
    % Calculate the adjusted Rand index
    ARI = (nij - ni*nj/n) / ((ni+nj)/2 - ni*nj/n);
end

function [Contingency, chi2, pval] = contingency(T1, T2)
    %CONTINGENCY - creates a contingency table for two vectors
    %
    % [Contingency, chi2,pval] = contingency(T1, T2)
    % This function takes two vectors of group assignments and creates a
    % contingency table. Additionally, it can perform a Chi-squared test
    % to test for independence between the clustering assignments.
    
    classes = unique([T1; T2]);
    Contingency = zeros(length(classes));
    
    for i = 1:length(classes)
        for j = 1:length(classes)
            Contingency(i,j) = sum(T1 == classes(i) & T2 == classes(j));
        end
    end
    
    % Perform chi-squared test
    [chi2, pval] = chi2test(Contingency);
end

function [chi2, pval] = chi2test(Contingency)
    %CHI2TEST - performs a Chi-squared test on a contingency table
    %
    % [chi2, pval] = chi2test(Contingency)
    % This function calculates the Chi-squared statistic and p-value for a
    % given contingency table.
    
    n = sum(Contingency(:)); % Total number of observations
    rowSums = sum(Contingency, 2);
    colSums = sum(Contingency, 1);
    
    expected = rowSums * colSums / n;
    
    chi2 = sum(sum((Contingency - expected).^2 ./ expected));
    df = (length(rowSums) - 1) * (length(colSums) - 1); % Degrees of freedom
    pval = 1 - chi2cdf(chi2, df);
end

