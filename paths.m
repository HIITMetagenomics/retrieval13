switch experiment
    case 'metaHIT',
        abundancePath = '/triton/ics/scratch/mi/ERA000116/abundance.mat';
        kmerPath = '/triton/ics/scratch/mi/ERA000116/kmerFrequency/ERA000116.k';
        kmerIndex = '';
        stringPath = '/triton/ics/scratch/mi/ERA000116/string.normalizedEntropy.C2/';
        stringPath = '/triton/ics/scratch/mi/ERA000116/string.between12n30.C2/';
        figfamPath = '/triton/ics/scratch/mi/ERA000116/figfam.normalizedEntropy.C2/';
        stringIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 25;
    case 'synthHigh',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e6r/abundance.mat';
        kmerPath = ['/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e6r/kmerFrequency/synth.high.k'];
        kmerIndex = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/multidist/rowIndex.txt';
        stringPath = ['/triton/ics/scratch/mi/metaSim/string-pHWFA-8g-46s-5e6r.normalizedEntropy.C2/'];
        figfamPath = ['/triton/ics/scratch/mi/metaSim/figfam-pHWFA-8g-46s-5e6r.normalizedEntropy.C2/'];
        stringIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 40;
    case 'synthLow',
        abundancePath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/abundance.mat';
        kmerPath = ['/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/synth.low.k'];
        kmerIndex = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/multidist/rowIndex.txt';
        stringPath = ['/triton/ics/scratch/mi/metaSim/string-pHWFA-8g-46s-5e5r.normalizedEntropy.C2/'];
        figfamPath = ['/triton/ics/scratch/mi/metaSim/figfam-pHWFA-8g-46s-5e5r.normalizedEntropy.C2/'];
        stringIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 57;
    case 'T2D-P1',
        abundancePath = '/triton/ics/scratch/mi/SRP008047/abundance.mat';
        clear label; load('/triton/ics/scratch/mi/SRP008047/annot_008047_matlab'); label = annot_008047_matlab(:,2); rowIndex = annot_008047_matlab(:,1); clear annot_008047_matlab
        kmerPath = '/triton/ics/scratch/mi/SRP008047/kmerFrequency/SRP008047.k';
        kmerIndex = '';
        stringPath = ['/triton/ics/scratch/mi/SRP008047/string.normalizedEntropy.C2/'];
        stringPath = ['/triton/ics/scratch/mi/SRP008047/string.between12n30.C2/'];
        figfamPath = ['/triton/ics/scratch/mi/SRP008047/figfam.normalizedEntropy.C2/'];
        stringIndex = '';
        entCutList = [0.01:0.01:1];
        posSample = 71;
    case 'T2D-P2'
        abundancePath = '/triton/ics/scratch/mi/SRP011011/abundance.mat';
        clear label; load('/triton/ics/scratch/mi/SRP011011/annot_011011_matlab'); label = annot_011011_matlab(:,2); rowIndex = annot_011011_matlab(:,1); clear annot_011011_matlab
        kmerPath = '/triton/ics/scratch/mi/SRP011011/kmerFrequency/SRP011011.k';
        kmerIndex = '/triton/ics/scratch/mi/SRP011011/rowIndex.txt';
        stringPath = ['/triton/ics/scratch/mi/SRP011011/string.normalizedEntropy.C2/'];
        stringPath = ['/triton/ics/scratch/mi/SRP011011/string.between12n30.C2/'];
        figfamPath = ['/triton/ics/scratch/mi/SRP011011/figfam.normalizedEntropy.C2/'];
        stringIndex = '';
        entCutList = [0.0:0.01:1];
        posSample = 99;
end

switch exp
    case 1
        subplot(2,2,3,'position',[0.1,0.1,0.35,0.35]);
    case 2
        subplot(2,2,1,'position',[0.1,0.6,0.35,0.35]);
    case 3
        subplot(2,2,2,'position',[0.6,0.6,0.35,0.35]);
    case 5
        subplot(2,2,4,'position',[0.6,0.1,0.35,0.35]);
end