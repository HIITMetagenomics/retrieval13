function changeWidthHorizontalLineErrorbar(figHandle,barLength)
% This function changes the width of the width of the errorbar
Xdata = get(figHandle,'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
% endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
%Xdata(xleft) = reshape((get(figHandle,'Xdata') * ones(1,length(Xdata(xleft))/length(get(figHandle,'Xdata'))))',length(Xdata(xleft)),1)' - barLength; %Xdata(xleft) - barLength;
%Xdata(xright) = reshape((get(figHandle,'Xdata') * ones(1,length(Xdata(xright))/length(get(figHandle,'Xdata'))))',length(Xdata(xright)),1)' + barLength; % Xdata(xright) + barLength;
Xdata(xleft) = reshape((ones(1,length(xleft)/length(get(figHandle,'Xdata')))' * get(figHandle,'Xdata')),length(Xdata(xleft)),1)' - barLength;
Xdata(xright) = reshape((ones(1,length(Xdata(xright))/length(get(figHandle,'Xdata')))' * get(figHandle,'Xdata')),length(Xdata(xright)),1)' + barLength;
set(figHandle,'Xdata',Xdata)
