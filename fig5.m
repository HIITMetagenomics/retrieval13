% Script for plotting figure
% Copyright Sohan Seth

clear all
load loadColorScheme; lists;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 6, 3]); POS = 0.5;

COMPUTE = 0;
if exist('figure5TestValues.mat','file')
    load figure5TestValues randomTest
else
    COMPUTE = 1;
end

for exp = [1,2,3,5]
    labelList = {};
    experiment = experimentList{exp}; paths;
    
    switch experiment
        case 'metaHIT',
            AXIS = [0 4 0 0.65];
        case 'synthHigh',
            AXIS = [0 4 0 0.95];
        case 'synthLow',
            AXIS = [0 4 0 0.9];
        case 'T2D-P1',
            AXIS = [0 4 0 0.7];
        case 'T2D-P2'
            AXIS = [0 4 0 0.75];
    end
    axis(AXIS)
    
    load('abundanceLabel.mat'); 
    abundanceLevelList = abundanceLevelListStore{exp};
    label = labelStore{exp};
    load([experiment,'.mat'])
    
    countLabel = 0; noSamples = length(label); posSample = sum(label);
    boxPlotData = zeros(posSample,4);
    d = zeros(noSamples,noSamples);
    avep = [];
    for count = 1:100
        avep = [avep; avgprec(d, label(randperm(length(label))))];
    end
    hold on, line([0, AXIS(2)], [mean(avep), mean(avep)],'linewidth',LINE,'color',randomColor,'linestyle','-'), hold off,
    line([0, AXIS(2)], [quantile(avep,QUP), quantile(avep,QUP)],'linewidth',LINE,'color',randomColor,'linestyle','--'), hold off,
    line([0, AXIS(2)], [quantile(avep,QDOWN), quantile(avep,QDOWN)],'linewidth',LINE,'color',randomColor,'linestyle','--'), hold off,
    aucBase = avep;

    for metric = {'count'}
        for M = 2
            MEAN = mean(squeeze(recordS(strcmp(metric, metricList), :, M == MList, :))');
            [maxVal, maxPos] = max(MEAN);
            labelList{end+1} = 'All';
            countLabel = countLabel +1; boxPlotData(:,countLabel) = squeeze(recordS(strcmp(metric, metricList), maxPos, M == MList, :));
        end
        
        if any(methodList == 2)
            for countK = kmerList
                for M = MList
                    MEAN = mean(squeeze(recordK(strcmp(metric, metricList), countK == kmerList, :, M == MList, :))');
                    [maxVal, maxPos] = max(MEAN);
                    hold off,
                    labelList{end+1} = num2str(countK);
                    countLabel = countLabel +1; boxPlotData(:,countLabel) = squeeze(recordK(strcmp(metric, metricList), countK == kmerList, maxPos, M == MList, :));
                end
            end
        end
    end
    
    if COMPUTE
        randomTest{exp} = zeros(4);
        for ii1 = 1:4
            for ii2 = ii1+1:4
                randomTest{exp}(ii1, ii2) = randomtest(boxPlotData(:, ii1), boxPlotData(:, ii2));
            end
        end
    end
    
    figure; hFig = distributionPlot(boxPlotData,'XValues',POS+(0:3),'showMM',2,'histOpt',1);
    options = struct('markersize',MARKERSIZE,'linewidth',LINE,'figHandle',hMain);
    distributionPlotHack(hFig,[stringColorList(1,:); kmerColorList],options);
    
    XX = get(get(gca,'child'),'Xdata'); XREF = XX{4}; XX = cell2mat(XX(3:-1:1));
    YY = get(get(gca,'child'),'Ydata'); YREF = YY{4}; YY = cell2mat(YY(3:-1:1));
    POSITION = get(gca,'position');
    for countArrow = 1:3
        if randomTest{exp}(1,countArrow+1) < 0.05
            if YY(countArrow) > YREF    
                text(POS+countArrow-0.04,AXIS(4),...
                    ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countArrow+1))))],...
                    'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom')                
            else
                text(POS+countArrow-0.04,AXIS(4),...
                    ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(1,countArrow+1))))],...
                    'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom')
            end
        end
    end
    set(gca,'xtick',(0:countLabel-1)+POS,'xticklabel',labelList,'fontsize',FONT);
    xlabel('$k$','fontsize',FONT,'interpreter','Latex')
    ylabel('Average precision','fontsize',FONT); box on
end

subplottitle
export_fig('../fig5.pdf');
if COMPUTE
    save figure5TestValues randomTest
end