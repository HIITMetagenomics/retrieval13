% Script for plotting figure
% Copyright Sohan Seth

clear all
load loadColorScheme; lists;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 6, 3]);  MARKERSIZE = 0.1;
POS = [0.3, 0.6]; 

COMPUTE = false;
if exist('figure6TestValues.mat','file')
    load figure6TestValues randomTest
else
    COMPUTE = true;
end

for exp = [1,2,3,5]
    pvalueCount = 1; 
    labelList = {};
    experiment = experimentList{exp}; paths;

    switch experiment
        case 'metaHIT',
            AXIS = [0 4 0.0 0.4];
        case 'synthHigh',
            AXIS = [0 4 0.0 0.7];
        case 'synthLow',
            AXIS = [0 4 0.0 0.65];
        case 'T2D-P1',
            AXIS = [0 4 0.48 0.7];
        case 'T2D-P2'
            AXIS = [0 4 0.0 0.7];
    end
    axis(AXIS)
    
    load('abundanceLabel.mat'); 
    abundanceLevelList = abundanceLevelListStore{exp};
    label = labelStore{exp};
    load([experiment,'.mat'])

    noSamples = length(label); posSample = sum(label);
    fprintf('Creating baseline from random distance matrices\n')
    d = zeros(noSamples,noSamples);
    avep = [];
    for count = 1:100
        avep = [avep; avgprec(d, label(randperm(length(label))))];
    end
    hold on, line([0, AXIS(2)], [mean(avep), mean(avep)],'linewidth',LINE,'color',randomColor,'linestyle','-'), hold off,
    
    countLabel = 0;
    for metric = {'count'}
         for M = MList
            MEAN = mean(squeeze(recordS(strcmp(metric, metricList), :, M == MList, :))');
            [maxVal, maxPos] = max(MEAN);
            avep = squeeze(recordS(strcmp(metric, metricList), [length(entCutList), maxPos], M == MList, :))';
            MEAN = mean(avep);
            STDE = std(avep)/sqrt(posSample);
            hold on, h = errorbar(countLabel + POS, MEAN, 0.5*STDE, 0.5*STDE, 'linewidth',LINE,'color',stringColorList(1,:),...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',stringColorList(1,:),'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.11);
            hold off,
            
            labelList{end+1} = ['All'];
            if COMPUTE
                [pvalue] = randomtest(avep(:,2), avep(:,1));
                randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
            end
            countLabel = countLabel +1;
            if randomTest{exp}(1,countLabel) < 0.05
                text(countLabel-1+POS(2)-0.04,AXIS(4),...
                    ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countLabel))))],...
                    'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',stringColorList(1,:))
            end
         end
        
        for countK = [12, 21, 30]
            for M = 2
                MEAN = mean(squeeze(recordK(strcmp(metric, metricList), countK == kmerList, :, M == MList, :))');
                [maxVal, maxPos] = max(MEAN);
                avep = squeeze(recordK(strcmp(metric, metricList), countK == kmerList, [length(entCutList), maxPos], M == MList, :))';
                MEAN = mean(avep);
                STDE = std(avep)/sqrt(posSample);
                hold on, h = errorbar(countLabel + POS, MEAN, 0.5*STDE, 0.5*STDE, 'linewidth',LINE,'color',kmerColorList(countK == kmerList,:),...
                    'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',kmerColorList(countK == kmerList,:),'markeredgecolor',[1 1 1]);
                changeWidthHorizontalLineErrorbar(h, 0.1);
                hold off,
                
                labelList{end+1} = num2str(countK);
                if COMPUTE
                    pvalue = randomtest(avep(:,2), avep(:,1));
                    randomTest{exp}(pvalueCount) = pvalue; pvalueCount = pvalueCount + 1;
                end
                countLabel = countLabel +1;
                if randomTest{exp}(1,countLabel) < 0.05
                    text(countLabel-1+POS(2)-0.04,AXIS(4),...
                        ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countLabel))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',kmerColorList(countK == kmerList,:))
                end
            end
        end
    end
    
    set(gca,'xtick',mean(POS) + (0:countLabel-1),'xticklabel',labelList,'fontsize',FONT);
    xlabel('$k$','fontsize',FONT,'interpreter','Latex');
    ylabel('Mean average precision','fontsize',FONT); box on   
end

subplottitle
export_fig(['../fig6.pdf']);
if COMPUTE
    save figure6TestValues randomTest
end