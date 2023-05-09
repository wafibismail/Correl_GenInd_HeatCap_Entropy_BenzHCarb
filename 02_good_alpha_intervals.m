% In this script, six plots are generated which shows α intervals for good
%  - correlation coefficient ρ between E and R_α
%  - correlation coefficient ρ between E and SCI_α
%  - correlation coefficient ρ between E and SO_α
%  - correlation coefficient ρ between ΔH and R_α
%  - correlation coefficient ρ between ΔH and SCI_α
%  - correlation coefficient ρ between ΔH and SO_α

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 3;
fontSize = 26;
saveToFile = false;
% Note: Dimensions of resulting images are scaled according to each window size.
%       To control the dimensions, after running the whole script, resize each
%       ... figure window and then run only the saveas functions
%       ... manually, by selection

% Cell containing Entropy and Heat Capacity of lower benzenoid
expData = {reshape([ % Entropy
    269.72 334.15 389.48 395.88 444.72 447.44 457.96 455.84 450.42 399.49
    499.83 513.86 508.54 507.39 506.08 512.52 500.73 510.31 509.21 513.88
    511.77 509.61 461.55 463.74 468.71 555.41 472.30 554.78 468.80 551.71
  ]', 30, 1), reshape([ % Heat Capacity
    83.019 133.325 184.194 183.654 235.165 233.497
   234.568 234.638 233.558 200.815 286.182 285.056
   284.037 284.088 285.148 284.595 284.870 284.503
   284.785 284.740 284.233 284.552 251.175 250.568
   251.973 336.098 267.543 337.204 285.041 368.518
  ]', 30, 1);
  "E", "ΔH" % Their labels
};

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 6 7  6 8  7 8 9 6  6  7  8  8  7 10 9  9  9  8  9  8 8 8  7 10  7  6  6  6
  0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8 10 10 10 12 10 12 8 8 10 12 10 20 12 16
  0 1 2 3  3 5  4 5 6 5  4  5  6  6  5 8  7  7  7  6  7  6 8 8  7  9 10  5 12 19
]'; % Used for computing indices based on edge endpoint degree partitions

% Cell containing the three index computing functions
getIndexFns = { % gets indices of all benzenoids, accepting variable alpha as argument
  @(a) (d_f(:,1)*4^a + d_f(:,2)*6^a  + d_f(:,3)*9^a); % General Randic index
  @(a) (d_f(:,1)*4^a + d_f(:,2)*5^a  + d_f(:,3)*6^a); % General SCI
  @(a) (d_f(:,1)*8^a + d_f(:,2)*13^a + d_f(:,3)*18^a) % General Sombor index
}';
% Cell containing their respective labels
indexName = {"R" "SCI" "SO"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two
numIndices = size(getIndexFns,2); % three
numCurves = numData*numIndices;   % six

% Boundaries for visible intervals, for each index-property pair
%             R_a   SCI_a   SO_a
xstart = [   -3.2    -5.5   -2.2   % E
             -2.4    -4.4   -1.9]; % ΔH
xend = [     -0.5    -1.5   -1.2   % E
             0.25    0.25   -0.1]; % ΔH
ystart = [0.99301 0.98901  0.983;  % E
           1.0001  1.0001 0.9981]; % ΔH
yend = [     0.95    0.96  0.975;  % E
             0.97   0.975   0.98]; % ΔH

% Boundaries for good alpha intervals
%                    R_a            SCI_a              SO_a
a_lb = [-2.5110455626637 -4.4900727614681 -1.96418016830175    % E
        -1.8384125797782 -3.3914310005448 -1.55591746776294];  % ΔH
a_ub = [-1.3332936534986 -2.7640513282069 -1.47956927233366    % E
        -0.54993661705396 -1.1480473895238 -0.60370868891694]; % ΔH
a_goodrho = [0.98; 0.99];

% Colors (different shades of cyan and green)
colShaded = {[0.85, 1, 1]; [0.85, 1, 0.85]};
colIndicatorVert = {[0.2, 0.55, 0.55]; [0, 0.5, 0]};
colIndicatorHorz = {[0.35, 0.75, 0.75]; [0.45, 0.7, 0.45]};
colCurve = {[0, 0.75, 0.75]; [0, 0.5, 0]};

% Do the same procedure for each experimental data i.e., E and ΔH
for ii = 1:numData
  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    ccFn = @(alpha) corrcoef( % Get: Corrcoef between benzenoid prop and index
      getIndexFns{n}(alpha)(!isnan(expData{1,ii})),
      expData{1,ii}(!isnan(expData{1,ii}))
    )(1,2);

    this_fig = figure((ii-1)*numIndices+n);
    hold on;

    % Plot the actual curve, not including good alpha range
    % generate x values
    x = [linspace(xstart(ii,n),a_lb(ii,n),300) linspace(a_ub(ii,n),xend(ii,n),300)];
    % generate their corresponding y values
    y = arrayfun(ccFn, x);
    plot(x, y, '-', 'LineWidth', lineWidth);

    % Shade the area between indicator lines
    x1 = linspace(a_lb(ii,n),a_ub(ii,n),400);
    a_goodrho_curve = arrayfun(ccFn, x1);
    iLine1 = zeros(1,size(x1,2)) + ystart(ii,n);
    iLine2 = zeros(1,size(x1,2)) + yend(ii,n);
    x2 = [x1, fliplr(x1)];
    inBetween = [iLine1, fliplr(iLine2)];
    fill(x2, inBetween, colShaded{ii}, 'LineStyle', 'none');

    % Draw the indicator lines
    plot([a_lb(ii,n) a_lb(ii,n)], [ystart(ii,n) yend(ii,n)], '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    plot([a_ub(ii,n) a_ub(ii,n)], [ystart(ii,n) yend(ii,n)], '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    plot([xstart(ii,n) a_lb(ii,n)], [a_goodrho(ii) a_goodrho(ii)], '--k', 'LineWidth', lineWidth/1.75);
    plot([a_lb(ii,n) a_ub(ii,n)], [a_goodrho(ii) a_goodrho(ii)], '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorHorz{ii});

    % Label the indicator lines
    text(a_lb(ii,n), yend(ii,n), {'', sprintf("  α=−%.04f", abs(a_lb)(ii,n))}, 'VerticalAlignment', 'top', 'Rotation', 90);
    text(a_ub(ii,n), yend(ii,n), {'', sprintf("  α=−%.04f", abs(a_ub)(ii,n))}, 'VerticalAlignment', 'top', 'Rotation', 90);

    % Plot the curve within the indicator lines with a different color
    plot(x1, a_goodrho_curve, '-', 'LineWidth', lineWidth, 'Color', colCurve{ii});
    plot(x1, iLine2, '-', 'LineWidth', lineWidth, 'Color', colCurve{ii});

    % Label the plot
    title(sprintf('between %s and %s_a', expData{2,ii}, indexName{n}));
    xlabel('α');
    ylabel('ρ');
    drawnow;
    axis([xstart(ii,n) xend(ii,n) yend(ii,n) ystart(ii,n)]);

    % Replace all hypens with minuses
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    % Set all fontsizes to size specified early in the script
    set(findall(this_fig,'-property','FontSize'),'FontSize', fontSize);
    drawnow;

    hold off;

  end
end

if saveToFile
  saveas(figure(1), "02_good_a_intervals_E_R_a.png");
  saveas(figure(2), "02_good_a_intervals_E_SCI_a.png");
  saveas(figure(3), "02_good_a_intervals_E_SO_a.png");
  saveas(figure(4), "02_good_a_intervals_DH_R_a.png");
  saveas(figure(5), "02_good_a_intervals_DH_SCI_a.png");
  saveas(figure(6), "02_good_a_intervals_DH_SO_a.png");
end
