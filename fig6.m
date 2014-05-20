% Script for plotting figure
% Copyright Sohan Seth

clearvars -except REAL useMetric
load loadBioColor; %loadColorScheme; 
lists;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 7.3, 2]);  MARKERSIZE = 0.1;
POS = [0.1, 0.3, 0.5]; addLayer = [0, 0.07];

COMPUTE = false;
if exist('figure6TestValues.mat','file')
    load figure6TestValues randomTest
else
    COMPUTE = true;
end

%REAL = 0;
if REAL
    plotExps = [11,5,7];
else
    plotExps = [6,8,9,10];
end
for exp = plotExps
    pvalueCount = 1; 
    labelList = {};
    experiment = experimentList{exp}; paths;
    if exp == 11;
        useMetric = {'log'}; %$$%
    else
        useMetric = {'count'};
    end
    titleList{exp} = createName(experiment, char(useMetric{1}));    
    
    switch experiment
        case 'metaHIT',
            AXIS = [0 5 0.0 0.4];
        case 'synthHigh',
            AXIS = [0 5 0.0 0.7];
        case 'synthLow',
            AXIS = [0 5 0.0 0.65];
        case 'T2D-P1',
            AXIS = [0 5 0.48 0.7];
        case 'T2D-P2'
            AXIS = [0 5 0.0 0.7];
        case 'bioRev'
            AXIS = [0 5 0.0 1];
        case 'bioRev-2'
            AXIS = [0 5 0.0 1];
        case 'HMP'
            AXIS = [0 5 0.0 1];
    end
    axis(AXIS)
    
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

    noSamples = length(label); posSample = sum(label);
    fprintf('Creating baseline from random distance matrices\n')
    d = zeros(noSamples,noSamples);
    avep = avgprect(d, label, testing);
    hold on, line([0, AXIS(2)], [mean(avep), mean(avep)],'linewidth',LINE,'color',randomColor,'linestyle','-'), hold off,
    
    countLabel = 0;
    for metric = useMetric%{'count'}
         for M = MList
            %MEAN = mean(squeeze(recordS(strcmp(metric, metricList), :, M == MList, :))');
            %[maxVal, maxPos] = max(MEAN);
            avep = [squeeze(recordS(strcmp(metric, metricList), length(entCutList), M == MList, :)), ...
                squeeze(recordS_cross(strcmp(metric, metricList), 1, M == MList, :)), ...
                squeeze(recordS_avgOverEntCut(strcmp(metric, metricList), 1, M == MList, :))];
            MEAN = mean(avep);
            STDE = std(avep)/sqrt(posSample);
            hold on, h = errorbar(countLabel + POS, MEAN, 0.5*STDE, 0.5*STDE, 'linewidth',LINE,'color',stringColorList(1,:),...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',stringColorList(1,:),'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.11);
            hold off,
            
            labelList{end+1} = ['All'];
            if COMPUTE
                for countOther = 2:size(avep,2)
                    [pvalue] = randomtest(avep(:,countOther), avep(:,1));
                    randomTest{exp}(countOther-1, pvalueCount) = pvalue;
                end
                %[pvalue] = randomtest(avep(:,3), avep(:,1));
                %randomTest2{exp}(pvalueCount) = pvalue;
                pvalueCount = pvalueCount + 1;
            end
            countLabel = countLabel +1;
            for countOther = 2:size(avep,2)
                if randomTest{exp}(countOther-1, countLabel) < 0.05 && mean(avep(:,countOther)) > mean(avep(:,1))
                    text(countLabel-1+POS(2)-0.04,AXIS(4)*(1 + addLayer(countOther-1)),...
                        ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(countOther-1,countLabel))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',stringColorList(1,:))
                end
                if randomTest{exp}(countOther-1, countLabel) < 0.05 && mean(avep(:,countOther)) < mean(avep(:,1))
                    text(countLabel-1+POS(2)-0.04,AXIS(4)*(1 + addLayer(countOther-1)),...
                        ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(countOther-1,countLabel))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',stringColorList(1,:))
                end
            end
            %if randomTest2{exp}(1,countLabel) < 0.05
            %    text(countLabel-1+POS(3)-0.04,AXIS(4),...
            %        ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countLabel))))],...
            %        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',stringColorList(1,:))
            %end
         end
        
        for countK = [12, 21, 30]
            for M = 2
                %MEAN = mean(squeeze(recordK(strcmp(metric, metricList), countK == kmerList, :, M == MList, :))');
                %[maxVal, maxPos] = max(MEAN);
                avep = [squeeze(recordK(strcmp(metric, metricList), countK == kmerList, length(entCutList), M == MList, :)), ...
                    squeeze(recordK_cross(strcmp(metric, metricList),  countK == kmerList, 1, M == MList, :)), ...
                    squeeze(recordK_avgOverEntCut(strcmp(metric, metricList),  countK == kmerList, 1, M == MList, :))];
                MEAN = mean(avep);
                STDE = std(avep)/sqrt(posSample);
                hold on, h = errorbar(countLabel + POS, MEAN, 0.5*STDE, 0.5*STDE, 'linewidth',LINE,'color',kmerColorList(countK == kmerList,:),...
                    'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',kmerColorList(countK == kmerList,:),'markeredgecolor',[1 1 1]);
                changeWidthHorizontalLineErrorbar(h, 0.1);
                hold off,
                
                labelList{end+1} = num2str(countK);
                if COMPUTE
                    for countOther = 2:size(avep,2)
                        [pvalue] = randomtest(avep(:,countOther), avep(:,1));
                        randomTest{exp}(countOther-1, pvalueCount) = pvalue;
                    end
                    pvalueCount = pvalueCount + 1;
                end
                countLabel = countLabel +1;
                for countOther = 2:size(avep,2)
                    if randomTest{exp}(countOther-1, countLabel) < 0.05 && mean(avep(:,countOther)) > mean(avep(:,1))
                        text(countLabel-1+POS(2)-0.04,AXIS(4)*(1 + addLayer(countOther-1)),...
                            ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(countOther-1,countLabel))))],...
                            'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',kmerColorList(countK == kmerList,:))
                    end
                    if randomTest{exp}(countOther-1, countLabel) < 0.05 && mean(avep(:,countOther)) < mean(avep(:,1))
                        text(countLabel-1+POS(2)-0.04,AXIS(4)*(1 + addLayer(countOther-1)),...
                            ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(countOther-1,countLabel))))],...
                            'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',kmerColorList(countK == kmerList,:))
                    end
                end
            end
        end
        
        for fig = 21
            %MEAN = mean(squeeze(recordK(strcmp(metric, metricList), countK == kmerList, :, M == MList, :))');
            %[maxVal, maxPos] = max(MEAN);
            avep = [squeeze(recordF(strcmp(metric, metricList), length(entCutList), fig == figList, :)), ...
                squeeze(recordF_cross(strcmp(metric, metricList),  1, fig == figList, :)), ...
                squeeze(recordK_avgOverEntCut(strcmp(metric, metricList), 1, fig == figList, :))];
            MEAN = mean(avep);
            STDE = std(avep)/sqrt(posSample);
            hold on, h = errorbar(countLabel + POS, MEAN, 0.5*STDE, 0.5*STDE, 'linewidth',LINE,'color',kmerColorList(countK == kmerList,:),...
                'linestyle','-','marker','s','markersize',MARKERSIZE,'markerfacecolor',kmerColorList(countK == kmerList,:),'markeredgecolor',[1 1 1]);
            changeWidthHorizontalLineErrorbar(h, 0.1);
            hold off,
            
            labelList{end+1} = 'FIG';
            if COMPUTE
                for countOther = 2:size(avep,2)
                    [pvalue] = randomtest(avep(:,countOther), avep(:,1));
                    randomTest{exp}(countOther-1, pvalueCount) = pvalue;
                end
                pvalueCount = pvalueCount + 1;
            end
            countLabel = countLabel +1;
            for countOther = 2:size(avep,2)
                if randomTest{exp}(countOther-1, countLabel) < 0.05 && mean(avep(:,countOther)) > mean(avep(:,1))
                    text(countLabel-1+POS(2)-0.04,AXIS(4)*(1 + addLayer(countOther-1)),...
                        ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(countOther-1,countLabel))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',kmerColorList(countK == kmerList,:))
                end
                 if randomTest{exp}(countOther-1, countLabel) < 0.05 && mean(avep(:,countOther)) < mean(avep(:,1))
                    text(countLabel-1+POS(2)-0.04,AXIS(4)*(1 + addLayer(countOther-1)),...
                        ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(countOther-1,countLabel))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',kmerColorList(countK == kmerList,:))
                end
            end
        end
        
    end
    
    set(gca,'xtick',mean(POS) + (0:countLabel-1),'xticklabel',labelList,'fontsize',FONT);
    xlabel('$k$','fontsize',FONT,'interpreter','Latex');  box on 
    if exp == plotExps(1)
        ylabel('Mean average precision','fontsize',FONT);  
    end
end

subplottitle;
saveEPS = false;
if saveEPS
    if REAL
        saveas(gcf,'../bioinformatics/fig6.eps','epsc')
    else
        saveas(gcf,'../bioinformatics/fig6a.eps','epsc')
    end
end
if COMPUTE
    save figure6TestValues randomTest
end