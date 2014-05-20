% order of data in Matlab matrices
experimentList = {'metaHIT', 'synthHigh', 'synthLow', 'T2D-P1', 'T2D-P2','bioRev','HMP','bioRev-2','bioRev-3','bioRev-4','metaHIT'}; % two metaHIT for plotting simplicity, synthHigh, synthLow, T2D-P1 are obsolete
metricList = {'count', 'sqrt', 'log'}; 
methodList = [1, 2, 3, 4, 5]; % 2 - k-mer based computation, 3 - distributed string mining based, 5 - figfam based
expList = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; 
kmerList = [12, 21, 30]; % length of k-mers
figList = [21, 30]; % length of k-mers from figfam
MList = 2; % cut-off (obsolete)

% Parameters
MARKERSIZE = 5; FONT = 8; LINE = 0.5; 
QUP = 0.95; QDOWN = 0.05; % quantiles
SEED = 1000; % random seed