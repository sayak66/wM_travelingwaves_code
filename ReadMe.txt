----------------------------------------

Code Descriptions:

Traveling waves in the prefrontal cortex during working memory
Sayak Bhattacharya, Scott L. Brincat, Mikael Lundqvist, Earl K. Miller
Correspondence: ekmiller@mit.edu


----------------------------------------
Main Code:


1. sb_spiral_curl_master2 - inputs LFP data and calculates circular-circular correlations

2. sb_spiral_curl_count - uses output from sb_spiral_curl_master2 and extracts the needed coefficients for further analysis

3. sb_spiral_count_2D_euclidean - rotating vs planar waves from coefficients obtained and simulated data

4. sb_spikes - analyzing spike data


-----------------------------------------

curl_planar01_sp24r_sp42r_sp24_sp42_noiseavg_fix_5curls.mat - coefficients obtained for simulated 72 types of waves . 1-32 are different directed planar waves. 33-72 are rotating waves of different directions and wavelengths (see Supplemental Figure 1).