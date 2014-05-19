% Script for plotting figure
% Copyright Sohan Seth

clearvars -except REAL
load loadBioColor; %loadColorScheme;
lists;
% hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 6, 3]); MARKERSIZE = 0.1;
hMain = figure('color',[1 1 1],'papertype','a4','PaperPosition', [1, 1, 7.3, 2]); MARKERSIZE = 0.1;

%REAL = 1;
if REAL
    plotExps = [11,5,7];
else
    plotExps = [6,8,9,10];
end

for exp = plotExps
    
    labelList = {};
    experiment = experimentList{exp}; paths;
    if exp == 11;
        useMetric = {'log'}; %$$%
    else
        useMetric = {'count'};
    end
    load([experiment,'_cross.mat'])
    
    titleList{exp} = createName(experiment,char(useMetric{1}));
    clear stringCountList labelList, countLabel = 1;
    
    switch experiment
        case 'metaHIT',
            AXIS = [0 1 0 15];
        case 'bioRev-4',
            AXIS = [0.7 1 0 15];
        case 'bioRev',
            AXIS = [0.7 1 0 15];
        case 'bioRev-2',
            AXIS = [0.7 1 0 15];
        case 'bioRev-3',
            AXIS = [0.7 1 0 15];
        case 'T2D-P1',
            AXIS = [0 1 0 15];
        case 'T2D-P2'
            AXIS = [0 1 0 15];
        case 'HMP'
            AXIS = [0 1 0 15];
    end
    axis(AXIS)
    
    hold on
    plot([0 0], [0 0],'linestyle','--','color',kmerColorList(1,:))
    plot([0 0], [0 0],'linestyle','--','color',kmerColorList(2,:))
    plot([0 0], [0 0],'linestyle','--','color',kmerColorList(3,:))
    plot([0 0], [0 0],'color',stringColorList(1,:))
    plot([0 0], [0 0],'linestyle',':','color',stringColorList(1,:))
    
    for countK = [12, 21, 30]
        dirName = [kmerPath,num2str(countK),'.q30.M2/stringCount.mat'];
        load(dirName,'stringCount');
        stringCountList(countLabel, :) = stringCount(2:end);
        if exp == 1 | exp == 11 | exp == 5
            plot(entCutList(2:end), log10(stringCountList(countLabel, :)), 'linewidth', LINE, 'color', kmerColorList(countLabel,:),'linestyle','--')
        else
            plot(entCutList, log10(stringCountList(countLabel, :)), 'linewidth', LINE, 'color', kmerColorList(countLabel,:),'linestyle','--')
        end
        ind = recordK_best(strcmp(char(useMetric{1}), metricList), countK == kmerList);
        plot(entCutList(ind),log10(stringCountList(countLabel, ind)),'s','markersize',5,...
            'markerfacecolor',kmerColorList(countLabel,:),'markeredgecolor',kmerColorList(countLabel,:));
        labelList{countLabel} = num2str(countK); countLabel = countLabel + 1;
    end
    
    dirName = [stringPath,'stringCount.mat'];
    load(dirName,'stringCount');
    stringCountList(countLabel, :) = stringCount(2:end);
    if exp == 1 | exp == 11 | exp == 5
            plot(entCutList(2:end), log10(stringCountList(countLabel, :)), 'linewidth', LINE, 'color', stringColorList(1,:))
        else
         plot(entCutList, log10(stringCountList(countLabel, :)), 'linewidth', LINE, 'color', stringColorList(1,:))
    end
    ind = recordS_best(strcmp(char(useMetric{1}), metricList));
    plot(entCutList(ind),log10(stringCountList(countLabel, ind)),'s','markersize',5,...
        'markerfacecolor',stringColorList(1,:),'markeredgecolor',stringColorList(1,:));
    labelList{countLabel} = 'All'; countLabel = countLabel + 1;
    
    dirName = [figfamPath,'21/stringCount.mat'];
    load(dirName,'stringCount');
    stringCountList(countLabel, :) = stringCount(2:end);
    if exp == 1 | exp == 11 | exp == 5
            plot(entCutList(2:end), log10(stringCountList(countLabel, :)), 'linewidth', LINE, 'color', stringColorList(1,:),'linestyle',':')
        else
            plot(entCutList, log10(stringCountList(countLabel, :)), 'linewidth', LINE, 'color', stringColorList(1,:),'linestyle',':')
    end
    ind = recordF_best(strcmp(char(useMetric{1}), metricList), 1, 21 == figList);
    plot(entCutList(ind),log10(stringCountList(countLabel, ind)),'s','markersize',5,...
        'markerfacecolor',stringColorList(1,:),'markeredgecolor',stringColorList(1,:));
    labelList{countLabel} = 'F'; countLabel = countLabel + 1;
            
    hold off
    axis(AXIS)
    set(gca,'fontsize',FONT)
    xlabel('Entropy threshold','fontsize',FONT); box on
    if exp == plotExps(1)
        ylabel('log_{10}(# of k-mers)','fontsize',FONT);
    end
    if exp == plotExps(end)%8 | exp == 7
        legend(labelList,'location','northeast','orientation','horizontal');
    end
 end

subplottitle;
saveEPS = false;
if saveEPS
    if REAL
        saveas(gcf,'../bioinformatics/fig8.eps','epsc')
    else
        saveas(gcf,'../bioinformatics/fig8a.eps','epsc')
    end
end