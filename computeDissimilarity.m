function D = computeDissimilarity(data, distance, verbose)

% This function computes the nSamxnSam dissimilarity matrix D given
%   data:   either a nSamxnSam interaction matrix (e.g. intersection/union)
%           or nFeatxnSam feature representation matrix, and
%   distance:   type of distance function, e.g., setover, hel
%   D:      nSamxnSam disimilairty matrix, 0 in the diagonal

% Author: Sohan Seth, sohan.seth@hhiit.fi

if nargin == 2
    verbose = 0; 
end
nSam = size(data,2);

switch distance
    case 'setover' % Percentage overlap between strings present
        if (norm(data - data','fro')) % If not symetric assume that only half the matrix is given
            data = data + data' - diag(diag(data));
        end
        if data(1,2) > data(1,1) % Determine if data is union or intersection matrix
            if verbose, fprintf('Computing %dx%d set overlap matrix from union matrix.\n',nSam,nSam); end
            U = data;
            I = repmat(diag(I),1,length(I)) + repmat(diag(I)',length(I),1) - U;
            D = 1 - I./(U+1); % +1 added for regulaization and avoid divide by zero
        else
            if verbose, fprintf('Computing %dx%d set overlap matrix from intersection matrix.\n',nSam,nSam); end;
            I = data;
            I = I + eye(length(data)); % Regularization to avoid zero strings
            U = repmat(diag(I),1,length(I)) + repmat(diag(I)',length(I),1) - I;
        end
        D = 1 - I./U;
    case 'hel' % Hellinger distance between probability measure
        if all(data(:) >= 0) & (all(data(:) <= 1))
            if verbose, fprintf('Computing %dx%d hellinger distance matrix from distributions.\n',nSam,nSam); end;
        else
            if verbose, fprintf('Computing %dx%d hellinger distance matrix from counts.\n',nSam,nSam); end;
            data = data ./ repmat(sum(data),size(data,1),1);
        end
        D = 1 - sqrt(data') * sqrt(data);
    case 'bc' % Bray Curtis distance
        if verbose, fprintf('Computing %dx%d Bray Curtis dissimilarity matrix from distributions.\n',nSam,nSam); end;
        for countSam1 = 1:nSam
            for countSam2 = 1:nSam
                D(countSam1,countSam2) = sum(abs(data(:,countSam1) - data(:,countSam2))) ...
                    / sum(data(:,countSam1) + data(:,countSam2));
            end
        end
    case 'd2S'
        if verbose, fprintf('Computing %dx%d d^2_S distance matrix from counts.\n',nSam,nSam); end;
        data = bsxfun(@minus, data, sum(data)/size(data,1));
        for count1 = 1:nSam
           for count2 = 1:nSam
               d(count1, count2) = ...
                   sum(data(:, count1) .* data(:, count2) ./ sqrt(data(:, count1).^2 + data(:, count2).^2)) ...
                   / sqrt(sum(data(:, count1).^2 ./ sqrt(data(:, count1).^2 + data(:, count2).^2))) ...
                   / sqrt(sum(data(:, count2).^2 ./ sqrt(data(:, count1).^2 + data(:, count2).^2)));
           end
        end
        D = 0.5 * (1 - d);
    otherwise
        error('dissimilarity not defined')
end

for count = 1:nSam
    D(count,count) = 0;
end
if any(D(:)<0) | any(D(:)>1)
    error('Dissimilarity is outside [0,1]')
end
if any(isinf(D(:))) | any(isnan(D(:))) | any(~isreal(D(:)))
    error('Inf, NaN, or imaginary values present')
end