% In this script, six values are closely approximated via golden section search
%  - correlation coefficient ρ between E and R_a with the best alpha value
%  - correlation coefficient ρ between E and SCI_a with the best alpha value
%  - correlation coefficient ρ between E and SO_a with the best alpha value
%  - correlation coefficient ρ between ΔH and R_a with the best alpha value
%  - correlation coefficient ρ between ΔH and SCI_a with the best alpha value
%  - correlation coefficient ρ between ΔH and SO_a with the best alpha value
% In addition, these curves are plotted

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2;
fontSize = 13;
saveToFile = false;

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

% Do the same procedure for each experimental data i.e., E and ΔH
for ii = 1:numData

  % Set visible ranges:

  % Upper and lower bounds for alpha (zoomed-in plot)
  xend = 0.25;
  if ii==1 xstart = -4; % between indices and E
  else xstart = -3; end % between indices and ΔH

  % Zoomed-in plot - Generate x values
  x = xstart:0.025:xend;
  yend = 0; % - to be assigned the highest final rho value among all indices

  % Upper and lower bounds for alpha (zoomed-out plot)
  xend_far = 40;
  xstart_far = -40;

  % Zoomed-out plot - Generate xfar values
  xfar = xstart_far:0.1:xend_far;

  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;

  % Indicate inverval for which R_a is the better indicator (before corr-curves)

  % WARNING: these xmeet1 values are hardcoded, computed separately
  xmeet1 = [
    -2.5787628478868;   % for E
    -1.6555622124036955 % for ΔH
  ](ii); % of these 2, use only one per iteration
  ymeet1 = [
    0.977567833002572; % for E
    0.994298524387525  % for ΔH
  ](ii);
  xmeet2 = 0;
  ymeet2 = [
    0.919815128725857; % for E
    0.974989505892947  % for ΔH
  ](ii);

  plot([xmeet1 xmeet1], [0 1], '--b', 'LineWidth', lineWidth);
  plot([0 xmeet1], [1 1], '--b', 'LineWidth', lineWidth);
  plot([0 0], [0 1], '--b', 'LineWidth', lineWidth);

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    ccFn = @(alpha) corrcoef( % Get: Corrcoef between benzenoid prop and index
      getIndexFns{n}(alpha)(!isnan(expData{1,ii})),
      expData{1,ii}(!isnan(expData{1,ii}))
    )(1,2);

    % generate their corresponding y values
    y = arrayfun(ccFn, x);

    % and their corresponding yfar values
    yfar = arrayfun(ccFn, xfar);

    % Compute peak values and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, xstart, xend, 1e-15))
    peakCorrCoeff = ccFn(peakAlpha)

    % Generate curve label
    curveLabels{n} = sprintf("%s and %s_a", expData{2,ii}, indexName{n});

    figure(ii); % One zoomed-out plot for each expData
    hold on;
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    curves(n) = plot(x, y, '-', 'LineWidth', lineWidth);

    % Show in plot where the peak is, and draw indicator lines
    plot([peakAlpha peakAlpha], [0 peakCorrCoeff], '--k', 'LineWidth', lineWidth/2);
    plot([xstart peakAlpha], [peakCorrCoeff peakCorrCoeff], '--k', 'LineWidth', lineWidth/2);
    text(peakAlpha, peakCorrCoeff,{'', sprintf("(−%.04f, %.04f)", abs(peakAlpha), peakCorrCoeff)}, 'VerticalAlignment', 'bottom');

    yend = max(yend, y(end)); % y value to be used as y lower bound
  end

  % Continue drawing indicator for interval for which R_a is best (text)
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b', 'MarkerSize', 16, 'LineWidth', lineWidth/1.5); % Mark them with an asterisk
  meet1Text = text(xmeet1, ymeet1, {'', sprintf(" (−%.04f, %.04f)", abs(xmeet1), ymeet1)}, 'VerticalAlignment', 'top');
  text(xmeet2, ymeet2, {'', sprintf("(0, %.04f) ", ymeet2)}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);
  set(meet1Text, 'Color', [0, 0, 0.8]); % Set to blue color, same as lines

  % Label this expData's zoomed-in plot
  xlabel('α');
  ylabel('ρ');
  leg = legend(curves, curveLabels);
  set(leg, 'location', "southwest");
  axis([xstart xend yend 1.001]);
  drawnow;

  % Label the zoomed-out plot
  figure(ii);
  xlabel('α');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels);
  if ii==2
    set(leg, 'location', "southeast");
  end

  hold off;
end

for ii = 1:4
  % Replace hyphen-minuses with hyphens
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));
end

for ii = 1:4
  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_E_indices_FAR.png");
  saveas(figure(2), "01_comb_ccurves_DH_indices_FAR.png");
  saveas(figure(3), "01_comb_ccurves_E_indices_NEAR.png");
  saveas(figure(4), "01_comb_ccurves_DH_indices_NEAR.png");
end
