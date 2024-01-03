% Number of false rejections divided by the total number of rejections.
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-10-17
% ------------------------------------------------------------------------------

function [isSig,sigLevel,pvals_corr] = BF_FDR(pvals,alpha,shutUp)
    
% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(alpha)
    % Default false discovery rate:
    alpha = 0.01;
end
if nargin < 3 || isempty(shutUp)
    shutUp = 1;
end

% ------------------------------------------------------------------------------
% Preliminaries:
% ------------------------------------------------------------------------------
if size(pvals,2) > size(pvals,1)
    pvals = pvals';
end

% Number of statistics tested, M
M = length(pvals);

% ------------------------------------------------------------------------------
% Sort p-values:
[pvals_sort, ix] = sort(pvals,'ascend');

% Probably a better way to do this but whatevs -- construct a reverse ordering:
iy = arrayfun(@(x)find(ix==x),1:M);

% Compute the number significant at the given alpha level:
numSig = find(pvals_sort < (1:M)'*alpha/M,1,'last');

% Compute the corrected p-values (alphas)
pvals_corr = M*pvals_sort./(1:M)';
pvals_corr = pvals_corr(iy); % sort backwards to original ordering

% Correct by averaging if there are duplicate p-values:
ups = unique(pvals);
if length(ups) < M
    for up = ups'
        dups = (pvals==up);
        if sum(dups)>1
            pvals_corr(dups) = mean(pvals_corr(dups))*ones(sum(dups),1);
        end
    end
end

% Return the significance level
if numSig > 0
    sigLevel = pvals_sort(numSig);
    if ~shutUp
        fprintf(1,'%u significant at p-threshold %.2g for a FDR=%.2g and %u samples\n',...
                                    numSig,sigLevel,alpha,M);
    end
else
    sigLevel = NaN;
    if ~shutUp
        fprintf(1,'***None significant for a FDR=%.2g and %u samples\n',alpha,M);
    end
end

% Return tests that are significant at the given alpha
isSig = logical(zeros(M,1));
isSig(ix(1:numSig)) = 1;

end