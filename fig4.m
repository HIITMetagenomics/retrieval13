% Script for plotting figure
% Copyright Sohan Seth

clearvars -except REAL useMetric
load loadBioColor; %loadColorScheme; 
lists;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 7.3, 2]); MARKERSIZE = 0.1;
POS = 0.5; 

COMPUTE = false;
if exist('figure4TestValues.mat','file')
    load figure4TestValues randomTest abundanceLevelListExp LOWCOV HIGHCOV kmerExp trueAbd
else
    COMPUTE = true;
end

if REAL
    plotExps = [1,5,7,11];
else
    plotExps = [6,8,9,10];
end
for exp = plotExps
    if exp == 11;
        useMetric = {'log'}; %$$%
    else
        useMetric = {'count'};
    end
    
    pvalueCount = 1; 
    labelList = {};
    experiment = experimentList{exp}; paths;
    titleList{exp} = createName(experiment, char(useMetric{1}));    
    
    switch experiment
        case 'metaHIT',
            AXIS = [0 4 0.0 0.4];
        case 'synthHigh',
            AXIS = [0 4 0.0 0.7];
        case 'bioRev',
            AXIS = [0 4 0.0 1];
        case 'bioRev-1',
            AXIS = [0 4 0.0 1];
        case 'bioRev-2',
            AXIS = [0 4 0.0 1];
        case 'bioRev-3',
            AXIS = [0 4 0.0 1];
        case 'synthLow',
            AXIS = [0 4 0.0 0.7];
        case 'T2D-P1',
            AXIS = [0 4 0.00 0.8];
        case 'T2D-P2'
            AXIS = [0 4 0.00 0.65];
        case 'HMP'
            AXIS = [0 4 0.00 1];
    end
    axis(AXIS)
    if any(exp == [2,3,6,8,9,10])
        axis(AXIS + [0 1 0 0])
    end
    
    if COMPUTE
        load(abundancePath)
        abundanceLevelListStore{exp} = bsxfun(@rdivide, abundanceLevelList, sum(abundanceLevelList,2));
        labelStore{exp} = label;
    else
        load('abundanceLabel.mat');
        abundanceLevelList = abundanceLevelListStore{exp};
        label = labelStore{exp};
    end
    load([experiment,'.mat'])
    load([experiment,'_avgOverEntCut.mat'])
    load([experiment,'_cross.mat'])
    
    randDiv = true;
    if randDiv
        rng(SEED); training = zeros(length(label), 1); testing = zeros(length(label), 1);
        temp = randperm(length(label)); training(temp(1:floor(end/2))) = 1; testing(temp(ceil(end/2):end)) = 1;
        training = logical(training);
    else
        testing = ones(length(label),1);
    end
    
    countLabel = 0;  posSample = sum(label); noSamples = length(label);
    d = zeros(noSamples,noSamples);
    avep = avgprect(d, label, testing);
    hold on, line([0, AXIS(2)+1], [mean(avep), mean(avep)],'linewidth',LINE,'color',randomColor,'linestyle','-'), hold off,
    aucBase = avep; % Baseline
    
    if any(exp == [2,3])
        TempAbd = abundanceLevelList;
        switch exp
            case 2
                abundanceLevelList = HIGHCOV';
            case 3
                abundanceLevelList = LOWCOV';
            otherwise
        end
    end
    
    %abundanceLevelList = abundanceLevelList ./ repmat(sum(abundanceLevelList,2),1,size(abundanceLevelList,2));
    %abundanceLevelListExp{exp} = abundanceLevelList;
    % abundanceLevelList = abundanceLevelListExp{exp}; %%% Updated to store
    abundanceLevelList = abundanceLevelListStore{exp};
    d = computeDissimilarity(abundanceLevelList', 'bc');
    aucB = avgprect(d, label, testing);
    figure(hMain); hold on, h = errorbar(AXIS(2)-POS -1, mean(aucB), ...
        0.5*sqrt(1/posSample)*std(aucB), 0.5*sqrt(1/posSample)*std(aucB), 'linewidth',LINE,'color',abundanceColor,...
        'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',abundanceColor,'markeredgecolor',[1 1 1]);
    changeWidthHorizontalLineErrorbar(h, 0.10);
    hold off,
    aucAbundance = aucB(:);
    
    if any(exp == [2,3])
        abundanceLevelList = TempAbd;
        abundanceLevelList = abundanceLevelList ./ repmat(sum(abundanceLevelList,2),1,size(abundanceLevelList,2));
        d = computeDissimilarity(abundanceLevelList', 'bc');
        aucB = avgprect(d, label, testing);
        figure(hMain); hold on, h = errorbar(AXIS(2)+POS, mean(aucB), ...
            0.5*sqrt(1/posSample)*std(aucB), 0.5*sqrt(1/posSample)*std(aucB), 'linewidth',LINE,'color',abundanceColor,...
            'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',abundanceColor,'markeredgecolor',[1 1 1]);
        changeWidthHorizontalLineErrorbar(h, 0.10);
        hold off,
        aucAbundance2 = aucB(:);
        abundanceLevelList = TempAbd;
    end
    
    if any(exp == [6,8,9,10])
        if COMPUTE
            TempAbd = abundanceLevelList;
            load(truePath,'abundanceLevelList');
            abundanceLevelList = abundanceLevelList ./ repmat(sum(abundanceLevelList,2),1,size(abundanceLevelList,2));
            d = computeDissimilarity(abundanceLevelList', 'bc');
            abundanceLevelList = TempAbd;
            trueAbd{exp} = d;
        else
            d = trueAbd{exp};
        end
        aucT = avgprect(d, label, testing);
        figure(hMain); hold on, h = errorbar(AXIS(2)+POS, mean(aucT), ...
            0.5*sqrt(1/posSample)*std(aucT), 0.5*sqrt(1/posSample)*std(aucT), 'linewidth',LINE,'color',abundanceColor,...
            'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',abundanceColor,'markeredgecolor',[1 1 1]);
        changeWidthHorizontalLineErrorbar(h, 0.10);
        hold off,
        aucAbundance2 = aucT(:);
    end

    if COMPUTE
        switch experiment
            case 'metaHIT',
                kmerPath = '/triton/ics/scratch/mi/ERA000116/kmerFrequency/ERA000116.k3.q30';
            case 'synthHigh',
                kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e6r/kmerFrequency/synth.high.k3.q30';
            case 'synthLow',
                kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/synth.low.k3.q30';
            case 'T2D-P1',
                kmerPath = '/triton/ics/scratch/mi/SRP008047/kmerFrequency/SRP008047.k3.q30';
            case 'T2D-P2'
                kmerPath = '/triton/ics/scratch/mi/SRP011011/kmerFrequency/SRP011011.k3.q30';
            case 'bioRev'
                kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r/kmerFrequency/bioRev.1e7r.k3.q30.append';
            case 'bioRev-2'
                kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-1e7r-2/kmerFrequency/bioRev-2.1e7r.k3.q30.append';
            case 'bioRev-3'
                kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-2e6r/kmerFrequency/bioRev-2e6r.k3.q30.append';
            case 'bioRev-4'
                kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-bioinformaticsRevision-mixed/kmerFrequency/bioRev-mixed.k3.q30.append';
            case 'HMP'
                kmerPath = '/triton/ics/scratch/mi/SRP002163/code/SRP002163.k3.q30';
        end
        if strcmp(experiment, 'HMP')
            [status, I] = system(['cut -d'' '' -f2- ',kmerPath]);
        else
            [status, I] = system(['cut -d'' '' -f2- ',kmerPath,' | grep -E -o ":[0-9]+" | grep -o -E "[0-9]+"']);
        end
        kmerExp{exp} = I;
    end
    I = kmerExp{exp};
    if strcmp(experiment, 'HMP')
        I = str2num(I); I = I';
    else
        I = str2num(I); I = reshape(I,length(I)/64,64);
    end
    %I = I ./ repmat(sum(I,2),1,size(I,2));
    d = computeDissimilarity(I', 'd2S');
    if strcmp(experiment,'T2D-P2')
        d_ = zeros(218); d_([1:158,160:197,199:218],[1:158,160:197,199:218]) = d;
        [~, I, ~] = intersect(load('rowIndexT2DP2'), rowIndex); d_ = d_(I, I); clear I; % kmerIndex
        aucB = avgprect(d_, label, testing);
    else
        aucB = avgprect(d, label, testing);
    end
    figure(hMain); hold on, h = errorbar(AXIS(2)-POS, mean(aucB), ...
        0.5*sqrt(1/posSample)*std(aucB), 0.5*sqrt(1/posSample)*std(aucB), 'linewidth',LINE,'color',k3Color,...
        'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',k3Color,'markeredgecolor',[1 1 1]);
    changeWidthHorizontalLineErrorbar(h, 0.10);
    hold off,
    aucK = aucB; % K-mer
    
    for metric = useMetric%{'log'}
        
        for fig = 21
            plotFig = 1;
            if plotFig
                % avep = squeeze(recordF_avgOverEntCut(strcmp(metric, metricList), 1, fig == figList, :));
                avep = squeeze(recordF(strcmp(metric, metricList), end, fig == figList, :));
                %avep = squeeze(recordF_cross(strcmp(metric, metricList), 1, fig == figList, :));
            else
                avep = aucBase;
            end
            MEAN = mean(avep);
            figure(hMain); hold on, h = errorbar(AXIS(2)-2-POS, MEAN, 0.5*sqrt(1/posSample)*std(avep), 0.5*sqrt(1/posSample)*std(avep), 'linewidth',LINE,'color',figColor,...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',figColor,'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.10);
            hold off,
            countLabel = countLabel +1;
        end
        aucFig = avep; % Figfam
        
        for M = 2
            %TEMP = squeeze(recordS_avgOverEntCut(strcmp(metric, metricList), 1, M == MList, :));
            TEMP = squeeze(recordS_cross(strcmp(metric, metricList), 1, M == MList, :));
            QUANUP = quantile(TEMP,QUP);
            QUANDOWN = quantile(TEMP,QDOWN);
            figure(hMain); hold on, h = errorbar(AXIS(2)-3-POS, mean(TEMP), 0.5*sqrt(1/posSample)*std(TEMP), ... 
                0.5*sqrt(1/posSample)*std(TEMP), 'linewidth',LINE,'color',stringColorList(1,:),...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',stringColorList(1,:),'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.10);
            hold off,
            avep = TEMP;
            countLabel = countLabel +1;
            if COMPUTE
                pvalue = randomtest(aucFig, avep);
                randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
                pvalue = randomtest(aucAbundance, avep);
                randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
                pvalue = randomtest(aucK, avep);
                randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
                if any(exp == [2,3])
                    pvalue = randomtest(aucAbundance2, avep);
                    randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
                end
                if any(exp == [6,8,9,10])
                    pvalue = randomtest(aucT, avep);
                    randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
                end
            end
        end
    end
    
    XX = get(get(gca,'child'),'Xdata'); XREF = XX{1}; XX = cell2mat(XX([2,4,3]));
    YY = get(get(gca,'child'),'Ydata'); YREF = YY{1}; YY = cell2mat(YY([2,4,3]));
    if any(exp == [2,3,6,8,9,10])
        XX = get(get(gca,'child'),'Xdata'); XREF = XX{1}; XX = cell2mat(XX([2,5,3,4]));
        YY = get(get(gca,'child'),'Ydata'); YREF = YY{1}; YY = cell2mat(YY([2,5,3,4]));
    end
    POSITION = get(gca,'position');
    for countArrow = 1:length(randomTest{exp})
        if randomTest{exp}(1,countArrow) < 0.05
            if YY(countArrow) > YREF
                text(countArrow+POS-0.04,AXIS(4),...
                    ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countArrow))))],...
                    'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom')
            else
                    text(countArrow+POS-0.04,AXIS(4),...
                        ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countArrow))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom')
            end
        end
    end
    
    labelList{1} = 'Ak'; 'All k';
    labelList{2} = 'FIG'; 'FIGfam';
    labelList{end+1} = 'Abd'; 'Abund.';
    labelList{end+1} = '3'; '3-mer';
    if any(exp == [2,3,6,8,9,10])
        labelList{end+1} = 'T'; 'True';
    end
    countLabel = length(labelList);
    set(gca,'xtick',(0:countLabel-1)+POS,'xticklabel',labelList,'fontsize',FONT)
    xlabel('Method','fontsize',FONT); box on
    if exp == plotExps(1)
        ylabel('Mean average precision','fontsize',FONT);
    end
end

subplottitle;
saveEPS = false;
if saveEPS
    if REAL
        saveas(gcf,'../bioinformatics/fig4.eps','epsc')
    else
        saveas(gcf,'../bioinformatics/fig4a.eps','epsc')
    end
end
if COMPUTE
    save figure4TestValues randomTest abundanceLevelListExp LOWCOV HIGHCOV kmerExp trueAbd
end