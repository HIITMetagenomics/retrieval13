% Script for plotting figure
% Copyright Sohan Seth

clear all
load loadColorScheme; lists;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 6, 3]); MARKERSIZE = 0.1;
POS = 0.5; 

COMPUTE = false;
if exist('figure4TestValues.mat','file')
    load figure4TestValues randomTest abundanceLevelListExp LOWCOV HIGHCOV kmerExp
else
    COMPUTE = true;
end

for exp = [1,2,3,5];
    pvalueCount = 1; 
    labelList = {};
    experiment = experimentList{exp}; paths;
    
    switch experiment
        case 'metaHIT',
            AXIS = [0 4 0.0 0.4];
        case 'synthHigh',
            AXIS = [0 4 0.0 0.7];
        case 'synthLow',
            AXIS = [0 4 0.0 0.7];
        case 'T2D-P1',
            AXIS = [0 4 0.00 0.8];
        case 'T2D-P2'
            AXIS = [0 4 0.00 0.65];
    end
    axis(AXIS)
    if any(exp == [2,3])
        axis(AXIS + [0 1 0 0])
    end
    
    load('abundanceLabel.mat'); 
    abundanceLevelList = abundanceLevelListStore{exp};
    label = labelStore{exp};
    load([experiment,'.mat'])
    
    countLabel = 0;  posSample = sum(label); noSamples = length(label);
    d = zeros(noSamples,noSamples);
    avep = [];
    for count = 1:100
        avep = [avep; avgprec(d, label(randperm(length(label))))];
    end
    hold on, line([0, AXIS(2)+1], [mean(avep), mean(avep)],'linewidth',LINE,'color',randomColor,'linestyle','-'), hold off,
    aucBase = avep;
    
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
    abundanceLevelList = abundanceLevelListExp{exp};
    d = computeDissimilarity(abundanceLevelList', 'hel');
    aucB = avgprect(d, label);
    figure(hMain); hold on, h = errorbar(AXIS(2)-POS -1, mean(aucB), ...
        0.5*sqrt(1/posSample)*std(aucB), 0.5*sqrt(1/posSample)*std(aucB), 'linewidth',LINE,'color',abundanceColor,...
        'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',abundanceColor,'markeredgecolor',[1 1 1]);
    changeWidthHorizontalLineErrorbar(h, 0.10);
    hold off,
    aucAbundance = aucB(:);
    
    if any(exp == [2,3])
        abundanceLevelList = TempAbd;
        abundanceLevelList = abundanceLevelList ./ repmat(sum(abundanceLevelList,2),1,size(abundanceLevelList,2));
        d = computeDissimilarity(abundanceLevelList', 'hel');
        aucB = avgprect(d, label);
        figure(hMain); hold on, h = errorbar(AXIS(2)+POS, mean(aucB), ...
            0.5*sqrt(1/posSample)*std(aucB), 0.5*sqrt(1/posSample)*std(aucB), 'linewidth',LINE,'color',abundanceColor,...
            'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',abundanceColor,'markeredgecolor',[1 1 1]);
        changeWidthHorizontalLineErrorbar(h, 0.10);
        hold off,
        aucAbundance2 = aucB(:);
        abundanceLevelList = TempAbd;
    end

%     switch experiment
%         case 'metaHIT',
%             kmerPath = '/triton/ics/scratch/mi/ERA000116/kmerFrequency/ERA000116.k3.q30';
%         case 'synthHigh',
%             kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e6r/kmerFrequency/synth.high.k3.q30';
%         case 'synthLow',
%             kmerPath = '/triton/ics/scratch/mi/metaSim/fasta/stringMining-pHWFA-8g-46s-5e5r/kmerFrequency/synth.low.k3.q30';
%         case 'T2D-P1',
%             kmerPath = '/triton/ics/scratch/mi/SRP008047/kmerFrequency/SRP008047.k3.q30';
%         case 'T2D-P2'
%             kmerPath = '/triton/ics/scratch/mi/SRP011011/kmerFrequency/SRP011011.k3.q30';
%     end
%     [status, I] = system(['cut -d'' '' -f2- ',kmerPath,' | grep -E -o ":[0-9]+" | grep -o -E "[0-9]+"']);
%     kmerExp{exp} = I;
    I = kmerExp{exp};
    I = str2num(I); I = reshape(I,length(I)/64,64);
    I = I ./ repmat(sum(I,2),1,size(I,2));
    d = computeDissimilarity(I', 'hel');
    if strcmp(experiment,'T2D-P2')
        d_ = zeros(218); d_([1:158,160:197,199:218],[1:158,160:197,199:218]) = d;
        [~, I, ~] = intersect(load('rowIndexT2DP2'), rowIndex); d_ = d_(I, I); clear I; % kmerIndex
        aucB = avgprect(d_, label);
    else
        aucB = avgprect(d, label);
    end
    figure(hMain); hold on, h = errorbar(AXIS(2)-POS, mean(aucB), ...
        0.5*sqrt(1/posSample)*std(aucB), 0.5*sqrt(1/posSample)*std(aucB), 'linewidth',LINE,'color',k3Color,...
        'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',k3Color,'markeredgecolor',[1 1 1]);
    changeWidthHorizontalLineErrorbar(h, 0.10);
    hold off,
    aucK = aucB;
    
    for metric = {'log'}
        for fig = 21
            MEAN = mean(squeeze(recordF(strcmp(metric, metricList), :, fig == figList, :))');
            [maxVal, maxPos] = max(MEAN);
            avep = squeeze(recordF(strcmp(metric, metricList), maxPos, fig == figList, :));
            MEAN = mean(avep);
            figure(hMain); hold on, h = errorbar(AXIS(2)-2-POS, MEAN, 0.5*sqrt(1/posSample)*std(avep), 0.5*sqrt(1/posSample)*std(avep), 'linewidth',LINE,'color',figColor,...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',figColor,'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.10);
            hold off,
            countLabel = countLabel +1;
        end
        aucFig = avep;
        
        for M = 2
            MEAN = mean(squeeze(recordS(strcmp(metric, metricList), :, M == MList, :))');
            [maxVal, maxPos] = max(MEAN);
            MEAN = MEAN([maxPos]);
            QUANUP = quantile(squeeze(recordS(strcmp(metric, metricList), [maxPos], M == MList, :))',QUP);
            QUANDOWN = quantile(squeeze(recordS(strcmp(metric, metricList), [maxPos], M == MList, :))',QDOWN);
            figure(hMain); hold on, h = errorbar(AXIS(2)-3-POS, MEAN, 0.5*sqrt(1/posSample)*std(squeeze(recordS(strcmp(metric, metricList), [maxPos], M == MList, :))'), 0.5*sqrt(1/posSample)*std(squeeze(recordS(strcmp(metric, metricList), [maxPos], M == MList, :))'), 'linewidth',LINE,'color',stringColorList(1,:),...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',stringColorList(1,:),'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.10);
            hold off,
            avep = squeeze(recordS(strcmp(metric, metricList), maxPos, M == MList, :));
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
            end
        end
    end
    
    XX = get(get(gca,'child'),'Xdata'); XREF = XX{1}; XX = cell2mat(XX([2,4,3]));
    YY = get(get(gca,'child'),'Ydata'); YREF = YY{1}; YY = cell2mat(YY([2,4,3]));
    if any(exp == [2,3])
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
    
    labelList{1} = 'All k';
    labelList{2} = 'FIGfam';
    labelList{end+1} = 'Abund.';
    labelList{end+1} = '3-mer';
    if any(exp == [2,3])
        labelList{end+1} = 'True';
    end
    countLabel = length(labelList);
    set(gca,'xtick',(0:countLabel-1)+POS,'xticklabel',labelList,'fontsize',FONT)
    xlabel('Method','fontsize',FONT)
    ylabel('Mean average precision','fontsize',FONT); box on
end

subplottitle;
export_fig('../fig4.pdf');
if COMPUTE
    save figure4TestValues randomTest abundanceLevelListExp LOWCOV HIGHCOV kmerExp labelExp
end