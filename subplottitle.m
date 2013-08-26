figure(hMain)
for countExp = [1,2,3,5]
    switch countExp
        case 1
            TITLE = ['(c) MetaHIT'];
            annotation('textbox',[0.00 0.50 0.2 0.01],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 2
            TITLE = ['(a) HIGH'];
            annotation('textbox',[0.00 0.99 0.2 0.01],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 3
            TITLE = ['(b) LOW'];
            annotation('textbox',[0.50 0.99 0.2 0.01],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
        case 5
            TITLE = ['(d) T2D-P2'];
            annotation('textbox',[0.50 0.50 0.2 0.01],'string',TITLE,'fontsize',FONT,'edgecolor',[1 1 1])
    end
end