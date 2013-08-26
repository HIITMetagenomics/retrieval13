function ap = avgprect(D, label)
% AVGPRECT returns the average precision values for each positive query
%   AP = PRECRECWRAP(D, LABEL) returns average precision values from 
%   distance matrix D for each positive sample as indicated in the 0/1 
%   LABEL vector. This function handles ties in scoring function using 
%   the method suggested by Frank McSherry and Marc Najork
%
%   Copyright Sohan Seth sohan.seth@hiit.fi

if length(D) ~= length(label)
    error('D and label dimension mismatch')
end

LOW = 10^6;
D = D - diag(LOW*ones(length(D),1)); % Force low diagonal entry
[~, sortInd] = sort(D(:,label == 1));
rel = label(sortInd(2:end,:)); % relevance binary values

[nSam, nQuery] = size(rel); % # of samples, # of queries
ap = zeros(sum(label), 1); % Average precision values
ind = (1:nSam); precision = zeros(size(ind));

for countQuery = 1:nQuery
    % Compute tc vector in the text, each element k is in (T(S(k)),T(S(k)+1)]
    D_ = D(:, countQuery); D_(D_ == -LOW) = [];
    [~, tc, S] = (unique(sort(D_)));
    tc = [0; tc];
    % Compute corresponding nc at each k
    nc = tc(2:end) - tc(1:end-1);
    rSample = sum(rel(:,countQuery)); % # of relevant samples for query
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
    ap(countQuery) = sum(precision(logical(rel(:,countQuery)))) / rSample;
end