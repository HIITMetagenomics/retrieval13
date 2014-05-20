function ap = avgprect(D, label, query)
% AVGPRECT returns the average precision values for each positive query
%   AP = AVGPRECT(D, LABEL) returns average precision values from 
%   distance matrix D for each positive sample as indicated in the 0/1 
%   LABEL vector. This function handles ties in scoring function using 
%   the method suggested by Frank McSherry and Marc Najork
%
%   UPDATE: If label is not binary but a set of integers then it queries with all samples
%   UPDATE: Only query with sample for which query is active
%
%   Copyright Sohan Seth sohan.seth@hiit.fi

if nargin == 2
    query = ones(size(label));
end

if size(D, 2) ~= length(label)
    error('D and label dimension mismatch')
end

HIGH = 10^6;
D = D - diag(HIGH*ones(length(D),1)); % Force HIGH diagonal entry

if min(size(label)) == 1 && max(unique(label(:))) == 1 % label is a binary vector, ...
    [~, sortInd] = sort(D(:,label == 1 & query == 1)); % ... query with only positive entries with query == true
    rel = label(sortInd(2:end,:)); % relevance binary values
    % ap = zeros(sum(label), 1); % Average precision values
else
    [~, sortInd] = sort(D(:, query == 1)); % ... query with all entries with query == true
    rel = bsxfun(@eq, label(sortInd(2:end,:)), label(logical(query))'); % relevance binary values
    % ap = zeros(length(label), 1); % Average precision values
end

[nSam, nQuery] = size(rel); % # of samples, # of queries
ind = (1:nSam); precision = zeros(size(ind));

ap = zeros(nQuery, 1);
for countQuery = 1:nQuery
    % Compute tc vector in the text, each element k is in (T(S(k)),T(S(k)+1)]
    D_ = D(:, countQuery); D_(D_ == -HIGH) = [];
    [~, tc, S] = (unique(sort(D_)));
    
    tc = [tc-1; length(D_)];
    
    % Compute corresponding nc at each k
    nc = tc(2:end) - tc(1:end-1);
    
    precisionTemp = cumsum(rel(:,countQuery));
    for countSam = 1:nSam
        if tc(S(countSam)) == 0
            precision(countSam) = countSam / nc(S(countSam)) ...
                * precisionTemp(tc(S(countSam)+1));
        else
        precision(countSam) = precisionTemp(tc(S(countSam))) + ...
            (countSam - tc(S(countSam))) / nc(S(countSam)) * ...
            (precisionTemp(tc(S(countSam)+1)) - precisionTemp(tc(S(countSam))));
        end
    end
    precision = precision(:) ./ ind(:); % precision @n
    
    % Average Precision, defined in section 2.4
    rSample = sum(rel(:,countQuery)); % # of relevant samples for query
    ap(countQuery) = sum(precision(logical(rel(:,countQuery)))) / rSample;
end
