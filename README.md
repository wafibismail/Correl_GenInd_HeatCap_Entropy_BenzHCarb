<h3>Files</h3>
<table>
  <tr><th>File Name</th><th>Description</th></tr>
  <tr><td>01_combined_corr_curves.m</td><td>Generates 4 combined plots for the correlation curves for the 6 pairs of data</td></tr>
  <tr><td>02_good_alpha_intervals.m</td><td>Generates 6 correlation curves, i.e., one for each pair of data, also indicating the alpha intervals which results in good correlation</td></tr>
  <tr><td>03_scatter_plots.m</td><td>Generates 6 scatter plots, one for each pair of data</td></tr>
  <tr><td>GoldenSectionSearch_Maximum.m</td><td>Octave 7.2 implementation of the <a href="https://en.wikipedia.org/wiki/Golden-section_search">Golden Section Search</a> algorithm. I manually translated this from the <a href="https://en.wikipedia.org/wiki/Golden-section_search">Python code</a>. This must be placed in the working directory as it is referenced in two of the other scripts.</td></tr>
</table>
<h3>Syntax</h3>
<p>The above scripts are to be run in Octave. Due to some syntax differences between Matlab and Octave, e.g. in subsetting matrices/arrays, the scripts will not run (without ammendments) in Matlab.</p>
<h3>Expected Output</h3>
<table>
  <tr><th colspan=6>Combined correlation curves</th></tr>
  <tr><td colspan=3><img src="imgs/01_comb_ccurves_DH_indices_FAR.png"></td><td colspan=3><img src="imgs/01_comb_ccurves_E_indices_FAR.png"></td></tr>
  <tr><td colspan=3><img src="imgs/01_comb_ccurves_DH_indices_NEAR.png"></td><td colspan=3><img src="imgs/01_comb_ccurves_E_indices_NEAR.png"></td></tr>
  <tr><th colspan=6>Good alpha intervals</th></tr>
  <tr><td><img src="imgs/02_good_a_intervals_DH_R_a.png"></td><td><img src="imgs/02_good_a_intervals_DH_SCI_a.png"></td><td><img src="imgs/02_good_a_intervals_DH_SO_a.png"></td><td><img src="imgs/02_good_a_intervals_E_R_a.png"></td><td><img src="imgs/02_good_a_intervals_E_SCI_a.png"></td><td><img src="imgs/02_good_a_intervals_E_SO_a.png"></td></tr>
  <tr><th colspan=6>Scatter plots</th></tr>
  <tr><td><img src="imgs/03_scatter_DH_R.png"></td><td><img src="imgs/03_scatter_DH_SCI.png"></td><td><img src="imgs/03_scatter_DH_SO.png"></td><td><img src="imgs/03_scatter_E_R.png"></td><td><img src="imgs/03_scatter_E_SCI.png"></td><td><img src="imgs/03_scatter_E_SO.png"></td></tr>
</table>
