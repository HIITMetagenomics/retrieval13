function ap = avgprec(D, label)
% AVGPREC returns the average precision values for each positive query
%   AP = PRECRECWRAP(D, LABEL) returns average precision values using
%   the 11 point interpolated precision recall from distance matrix D for
%   each positive sample as indicated in the 0/1 LABEL vector.
%
%   Copyright Sohan Seth sohan.seth@hiit.fi

D = D - diag(10^6*ones(length(D),1)); % Force low diagonal entry
if length(D) ~= length(label)
    error('D and label dimension mismatch')
end

[~, sortInd] = sort(D(:,label == 1));
rel = label(sortInd(2:end,:));
[nSam, nQuery] = size(rel); % # of samples, # of queries
ap = zeros(sum(label), 1); % Average precision values
iap11 = zeros(11,1); ind = (1:nSam); recPos = 0:0.1:1;
for countQuery = 1:nQuery
    rSample = sum(rel(:,countQuery)); % # of relevant samples for query
    precision = cumsum(rel(:,countQuery)) ./ ind(:); % precision @n
    recall = cumsum(rel(:,countQuery)) / rSample;
    for countRecall = recPos
        iap11(countRecall == recPos)  = max(precision(recall >= countRecall));
    end
    % Smoothing to reduce wiggles
    % ap(countQuery) = mean(iap11);
    % Direct computation without smoothing
    ap(countQuery) = sum(precision(logical(rel(:,countQuery)))) / rSample;
end