figure(hMain)
for countExp = plotExps
    switch countExp
        case 1
            TITLE = titleList{1};
            annotation('textbox',[0.13 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 9
            TITLE = titleList{9};
            annotation('textbox',[0.61 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 6
            TITLE = titleList{6};
            annotation('textbox',[0.13 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 8
            TITLE = titleList{8};
            annotation('textbox',[0.37 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 2
            TITLE = titleList{2};
            annotation('textbox',[0.13 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 7
            TITLE = titleList{7};
            annotation('textbox',[0.85 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 3
            TITLE = titleList{3};
            annotation('textbox',[0.37 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 5
            TITLE = titleList{5};
            annotation('textbox',[0.61 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 10
            TITLE = titleList{10};
            annotation('textbox',[0.85 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 11
            TITLE = titleList{11};
            annotation('textbox',[0.37 0.97 0.2 0.01] + [-0.03 0 0 0],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
    end
end