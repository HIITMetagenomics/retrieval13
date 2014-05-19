function name = createName(experiment, metric)

switch experiment
    case 'metaHIT'
        expName = 'MetaHIT';
    case 'HMP'
        expName = 'HMP';
    case 'T2D-P2'
        expName = 'T2D-P2';
    case 'bioRev'
        expName = 'HIGH-C';
    case 'bioRev-2'
        expName = 'HIGH-VAR';
    case 'bioRev-3'
        expName = 'LOW-C';
    case 'bioRev-4'
        expName = 'MIXED-C';
end

if nargin == 1
    name = [expName];
else
    if strcmp(metric, 'count')
        metricName = 'Jaccard';
    else
        metricName = metric;
    end
    name = [expName,' + ',metricName];
end