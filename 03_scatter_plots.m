% In this script, six scatter plots are generated for the following pairs
%  - E and R_{-1.9482}
%  - E and SCI_{-3.6624}
%  - E and SO_{-1.7285}
%  - ΔH and R_{-1.2383}
%  - ΔH and SCI_{-2.3554}
%  - ΔH and SO_{-1.1165}
% In addition, the curves of ρ against α and scatter plots of best indices are plotted

close all; % Before drawing, close any figures already opened
clear;     % Clear all variables

% CONFIG: Line width and font fize for each curve in drawn figures
lineWidth = 2;
fontSize = 26;
% Save plots to images? Set to true if yes
saveToFile = false;

% Cell containing Entropy and Heat Capacity of lower benzenoid
expData = {reshape([
    269.72 334.15 389.48 395.88 444.72 447.44 457.96 455.84 450.42 399.49
    499.83 513.86 508.54 507.39 506.08 512.52 500.73 510.31 509.21 513.88
    511.77 509.61 461.55 463.74 468.71 555.41 472.30 554.78 468.80 551.71
  ]', 30, 1), reshape([
    83.019 133.325 184.194 183.654 235.165 233.497
   234.568 234.638 233.558 200.815 286.182 285.056
   284.037 284.088 285.148 284.595 284.870 284.503
   284.785 284.740 284.233 284.552 251.175 250.568
   251.973 336.098 267.543 337.204 285.041 368.518
  ]', 30, 1);
  "E", "ΔH" % Their labels
}

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 6 7  6 8  7 8 9 6  6  7  8  8  7 10 9  9  9  8  9  8 8 8  7 10  7  6  6  6
  0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8 10 10 10 12 10 12 8 8 10 12 10 20 12 16
  0 1 2 3  3 5  4 5 6 5  4  5  6  6  5 8  7  7  7  6  7  6 8 8  7  9 10  5 12 19
]';

% Cell containing the three index computing functions
getIndexFns = { % gets indices of all benzenoids, accepting variable alpha as argument
  @(a) (d_f(:,1)*4^a + d_f(:,2)*6^a  + d_f(:,3)*9^a); % General Randic aka PCI
  @(a) (d_f(:,1)*4^a + d_f(:,2)*5^a  + d_f(:,3)*6^a); % General SCI
  @(a) (d_f(:,1)*8^a + d_f(:,2)*13^a + d_f(:,3)*18^a) % General Sombor index
}';
% Cell containing their respective labels
indexName = {"R" "SCI" "SO"};

numData = size(expData, 2);       % two
numIndices = size(getIndexFns,2); % three
numCurves = numData*numIndices;   % six

for edn = 1:numData
  for fnn = 1:numIndices
    ccFn = @(alpha) corrcoef( % Gets corrcoef between benzenoid prop and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);

    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));

    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);
    x = [getIndexFns{fnn}(peakAlpha)(1) getIndexFns{fnn}(peakAlpha)(end)];
    y = m*x + b;

    % Scatter plot
    this_figure = figure(3*(edn-1)+fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth', lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s_{−%.4f}", indexName{fnn}, abs(peakAlpha));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});

    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");

    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)

    drawnow;
  end
end

if saveToFile
  saveas(figure(1), "03_scatter_E_R.png");
  figure(2);
  axis([0.035 0.1125 250 600]);
  saveas(figure(2), "03_scatter_E_SCI.png");
  saveas(figure(3), "03_scatter_E_SO.png");
  saveas(figure(4), "03_scatter_DH_R.png");
  figure(5);
  axis([0.2 0.9625 50 400]);
  saveas(figure(5), "03_scatter_DH_SCI.png");
  saveas(figure(6), "03_scatter_DH_SO.png");
end
