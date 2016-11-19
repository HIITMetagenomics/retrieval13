function H = distributionPlotHack(figHandle,Colors,options)
% DISTRIBUTIONPLOTHACK plots the violin plot with smoothed density
%   boundaries rather than histograms as in DISTRIBUTIONPLOT
%   H = distributionPlotHack(figHandle,Colors,options) takes the handles
%   returned by distributionPlot as input, along with a set of Colors, and
%   optional arguments, and plot similar plot but with smoothed boundaries
%   and return the function handle.
%
% Author Sohan Seth sohan.seth@hiit.fi

if nargin < 3
    options.markersize = 6;
    options.linewdith = 2;
    options.figHandle = figure();
end
MARKERSIZE = options.markersize;
LINE = options.linewidth;
H = options.figHandle;

if nargin < 2
    Colors = [0 1 0] * ones(1, length(figHandle{1}));
end

figure(H)
hold on;
for count = 1:length(figHandle{1})
    X = get(figHandle{1}(count),'Xdata');
    Y = get(figHandle{1}(count),'Ydata');
    plot(X(1,:),Y(1,:),'color',Colors(count,:),'linewidth',LINE)
    plot(X(2,:),Y(2,:),'color',Colors(count,:),'linewidth',LINE)
end

X = get(figHandle{2},'Xdata');
Y = get(figHandle{2},'Ydata');
for count = 1:length(X)
    plot(X(count),Y(count),'s','markersize',MARKERSIZE,...
        'markerfacecolor',Colors(count,:),'markeredgecolor',Colors(count,:))
end
hold off
