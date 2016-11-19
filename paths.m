switch experiment
    case 'metaHIT',
        abundancePath = '/triton/ics/scratch/mi/ERA000116/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/ERA000116/kmerFrequency/ERA000116.k';
        kmerIndex = '';
        stringPath = '/triton/ics/scratch/mi/ERA000116/string.normalizedEntropy.C2/';
        stringPath = '/triton/ics/scratch/mi/ERA000116/string.between12n30.C2/';
        figfamPath = '/triton/ics/scratch/mi/ERA000116/figfam.normalizedEntropy.C2/';
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 25;
    case 'synthHigh',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e6r/abundance.mat';
        kmerPath = ['/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e6r/kmerFrequency/synth.high.k'];
        kmerIndex = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/multidist/rowIndex.txt';
        stringPath = ['/triton/ics/scratch/mi/metaSim/string-pHWFA-8g-46s-5e6r.normalizedEntropy.C2/'];
        figfamPath = ['/triton/ics/scratch/mi/metaSim/figfam-pHWFA-8g-46s-5e6r.normalizedEntropy.C2/'];
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 40;
    case 'synthLow',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/abundance.mat';
        kmerPath = ['/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/synth.low.k'];
        kmerIndex = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/multidist/rowIndex.txt';
        stringPath = ['/triton/ics/scratch/mi/metaSim/string-pHWFA-8g-46s-5e5r.normalizedEntropy.C2/'];
        figfamPath = ['/triton/ics/scratch/mi/metaSim/figfam-pHWFA-8g-46s-5e5r.normalizedEntropy.C2/'];
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 57;
    case 'T2D-P1',
        abundancePath = '/triton/ics/scratch/mi/SRP008047/abundance.mat';
        clear label; load('/triton/ics/scratch/mi/SRP008047/annot_008047_matlab'); label = annot_008047_matlab(:,2); rowIndex = annot_008047_matlab(:,1); clear annot_008047_matlab
        kmerPath = '/triton/ics/scratch/mi/SRP008047/kmerFrequency/SRP008047.k';
        kmerIndex = '';
        stringPath = '/triton/ics/scratch/mi/SRP008047/string.normalizedEntropy.C2/';
        stringPath = '/triton/ics/scratch/mi/SRP008047/string.between12n30.C2/';
        figfamPath = '/triton/ics/scratch/mi/SRP008047/figfam.normalizedEntropy.C2/';
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.01:0.01:1];
        posSample = 71;
    case 'T2D-P2'
        abundancePath = '/triton/ics/scratch/mi/SRP011011/abundance.mat';
        clear label; % load('/triton/ics/scratch/mi/SRP011011/annot_011011_matlab'); 
        load annot_011011_matlab; label = annot_011011_matlab(:,2); rowIndex = annot_011011_matlab(:,1); clear annot_011011_matlab
        kmerPath = '/triton/ics/scratch/mi/SRP011011/kmerFrequency/SRP011011.k';
        kmerIndex = '/triton/ics/scratch/mi/SRP011011/rowIndex.txt';
        stringPath = '/triton/ics/scratch/mi/SRP011011/string.normalizedEntropy.C2/';
        stringPath = '/triton/ics/scratch/mi/SRP011011/string.between12n30.C2/';
        figfamPath = '/triton/ics/scratch/mi/SRP011011/figfam.normalizedEntropy.C2/';
        figfamIndex = '';
        stringIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 99;
    case 'bioRev',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r/abundanceMetaPhlAn.mat';
        truePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r/kmerFrequency/bioRev-1e7r.k';
        kmerIndex = '';
        stringPath = ['/triton/ics/scratch/mi/metaSim/string-bioinformaticsRevision-1e7r.normalizedEntropy.C2/'];
        figfamPath = '/triton/ics/scratch/mi/metaSim/stringMining-bioinformaticsRevision-1er7.figfam/';
        stringIndex = '';
        figfamIndex = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r/stringMining-bioinformaticsRevision-1er7.figfam.txt';
        entCutList = [0.01:0.01:1];
        posSample = 98;
    case 'bioRev-2',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r-2/abundanceMetaPhlAn.mat';
        truePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r-2/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r-2/kmerFrequency/bioRev-1e7r.k';
        kmerIndex = '';
        stringPath = ['/triton/ics/scratch/mi/metaSim/string-bioinformaticsRevision-1e7r-2.normalizedEntropy.C2/'];
        figfamPath = '/triton/ics/scratch/mi/metaSim/stringMining-bioinformaticsRevision-1er7-2.figfam/';
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.01:0.01:1];
        posSample = 98;
    case 'bioRev-3',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-2e6r/abundanceMetaPhlAn.mat';
        truePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-2e6r/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-2e6r/kmerFrequency/bioRev-2e6r.k';
        kmerIndex = '';
        stringPath = '/triton/ics/scratch/mi/metaSim/string-bioinformaticsRevision-2e6r.normalizedEntropy.C2/';
        figfamPath = '/triton/ics/scratch/mi/metaSim/stringMining-bioinformaticsRevision-2e6r.figfam/';
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.01:0.01:1];
        posSample = 98;
    case 'bioRev-4',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-mixed/abundanceMetaPhlAn.mat';
        truePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-mixed/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-mixed/kmerFrequency/bioRev-mixed.k';
        kmerIndex = '';
        stringPath = '/triton/ics/scratch/mi/metaSim/string-bioinformaticsRevision-mixed.normalizedEntropy.C2/';
        figfamPath = '/triton/ics/scratch/mi/metaSim/stringMining-bioinformaticsRevision-mixed.figfam/';
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.01:0.01:1];
        posSample = 98;
    case 'HMP',
        abundancePath = '/triton/ics/scratch/mi/SRP002163/code/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/SRP002163/kmerFrequency/SRP002163.k';
        kmerIndex = '';
        stringPath = '/triton/ics/scratch/mi/SRP002163/string.normalizedEntropy.C2/';
        figfamPath = '/triton/ics/scratch/mi/SRP002163/figfam.normalizedEntropy.C2/';
        stringIndex = '';
        figfamIndex = '';
        entCutList = [0.01:0.01:1];
        posSample = 435;
end

% switch exp
%     case 1
%         subplot(2,2,3,'position',[0.1,0.1,0.35,0.35]);
%     case 2
%         subplot(2,2,1,'position',[0.1,0.6,0.35,0.35]);
%     case 3
%         subplot(2,2,2,'position',[0.6,0.6,0.35,0.35]);
%     case 5
%         subplot(2,2,4,'position',[0.6,0.1,0.35,0.35]);
% end

switch exp
    case 1
        subplot(1,4,3,'position',[0.10,0.3,0.165,0.5]);
    case 11
        subplot(1,4,2,'position',[0.34,0.3,0.165,0.5]);
    case 9
        subplot(1,4,3,'position',[0.58,0.3,0.165,0.5]);
    case 10
        subplot(1,4,4,'position',[0.82,0.3,0.165,0.5]);
    case 6
        subplot(1,4,1,'position',[0.10,0.3,0.165,0.5]);
    case 8
        subplot(1,4,2,'position',[0.34,0.3,0.165,0.5]);
    case 2
        subplot(1,4,1,'position',[0.10,0.3,0.165,0.5]);
    case 7
        subplot(1,4,4,'position',[0.82,0.3,0.165,0.5]);
    case 3
        subplot(1,4,2,'position',[0.34,0.3,0.165,0.5]);
    case 5
        subplot(1,4,3,'position',[0.58,0.3,0.165,0.5]);
end
