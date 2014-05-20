% Script for plotting figure
% Copyright Sohan Seth

clearvars -except REAL useMetric
load loadBioColor; %loadColorScheme; 
lists;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 7.3, 2]); POS = 0.5;

COMPUTE = false;
if exist('figure7TestValues.mat','file')
    load figure7TestValues randomTest
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
    labelList = {};
    experiment = experimentList{exp}; paths;
    
    switch experiment
        case 'metaHIT',
            AXIS = [0 3 0.00 0.65];
        case 'synthHigh',
            AXIS = [0 3 0.00 0.95];
        case 'synthLow',
            AXIS = [0 3 0.00 1.00];
        case 'T2D-P1',
            AXIS = [0 3 0.43 0.7];          
        case 'T2D-P2'
            AXIS = [0 3 0.0 0.75];
        case 'HMP'
            AXIS = [0 3 0.00 1];
        case 'bioRev',
            AXIS = [0 3 0.0 1];
        case 'bioRev-1',
            AXIS = [0 3 0.0 1];
        case 'bioRev-2',
            AXIS = [0 3 0.0 1];
        case 'bioRev-3',
            AXIS = [0 3 0.0 1];
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
    
    countLabel = 0; noSamples = length(label); posSamples = sum(label);
    clear boxPlotData % = zeros(posSample,length(metricList));
    
    d = zeros(noSamples,noSamples);
    avep = avgprect(d, label, testing);
    hold on, line([-1, AXIS(2)], [mean(avep), mean(avep)],'linewidth',LINE,'color',randomColor,'linestyle','-'), hold off,
    hold on, line([-1, AXIS(2)], [quantile(avep,QUP), quantile(avep,QUP)],'linewidth',LINE,'color',randomColor,'linestyle','--'), hold off,
    hold on, line([-1, AXIS(2)], [quantile(avep,QDOWN), quantile(avep,QDOWN)],'linewidth',LINE,'color',randomColor,'linestyle','--'), hold off,
    aucBase = avep;

    for metric = metricList
        titleList{exp} = createName(experiment);
        
        for M = 2;
            %MED = mean(squeeze(recordS(strcmp(metric, metricList), :, M == MList, :))');
            %[maxVal, maxPos] = max(MED);
            labelList{end+1} = 'All';
            avep = squeeze(recordS_cross(strcmp(metric, metricList), 1, M == MList, :)); %squeeze(recordS(strcmp(metric, metricList), maxPos, M == MList, :))';
            boxPlotData(:,strcmp(metric, metricList)) = avep(:);
            countLabel = countLabel +1;
        end
    end
    
    if COMPUTE
            randomTest{exp} = zeros(3);
            for ii1 = 1:3
                for ii2 = ii1+1:3
                    randomTest{exp}(ii1, ii2) = randomtest(boxPlotData(:, ii1), boxPlotData(:, ii2));
                end
            end
    end
    
    h = figure; hFig = distributionPlot(boxPlotData,'XValues',POS + (0:2),'showMM',2,'histOpt',1);
    options = struct('markersize',MARKERSIZE,'linewidth',LINE,'figHandle',hMain);
    distributionPlotHack(hFig,[kmerColorList;stringColorList(1,:)],options); close(h);
    
    arrowColors = [kmerColorList;stringColorList(1,:)];
    for i1 = 1:3
        i3 = 0;
        randomTest{exp} = (randomTest{exp} + randomTest{exp}')/2;
        for i2 = 1:3
            if i1 == i2 
                continue
            else
                i3 = i3 + 1;
            if randomTest{exp}(i1,i2) < 0.05
                if mean(boxPlotData(:,i1)) < mean(boxPlotData(:,i2))
                    text(POS-1+i1,AXIS(4)+(i3-1)*0.07,...
                        ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(i1,i2))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',arrowColors(i2,:));
                else
                    text(POS-1+i1,AXIS(4)+(i3-1)*0.07,...
                        ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(i1,i2))))],...
                        'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',arrowColors(i2,:));
                end
%                 if mean(boxPlotData(:,i1)) < mean(boxPlotData(:,i2))
%                     text(POS-1+i1-0.6+i3*0.3,AXIS(4),...
%                         ['\downarrow ',char(42*ones(1,decideStars(randomTest{exp}(i1,i2))))],...
%                         'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',arrowColors(i2,:));
%                 else
%                     text(POS-1+i1-0.6+i3*0.3,AXIS(4),...
%                         ['\uparrow ',char(42*ones(1,decideStars(randomTest{exp}(i1,i2))))],...
%                         'edgecolor',[1 1 1],'fontsize',FONT,'interpreter','tex','verticalalignment','bottom','color',arrowColors(i2,:));
%                 end
            end
            end
        end
    end
       
    figure(hMain)
    set(gca,'xtick',(0:2)+POS,'xticklabel',metricList,'fontsize',FONT)
    xlabel('Metric','fontsize',FONT); box on
    if exp == plotExps(1)
        ylabel('Average precision','fontsize',FONT);
    end
end

subplottitle
saveEPS = false;
if saveEPS
    if REAL
        saveas(gcf,'../bioinformatics/fig7.eps','epsc')
    else
        saveas(gcf,'../bioinformatics/fig7a.eps','epsc')
    end
end
if COMPUTE
    save figure7TestValues randomTest
end