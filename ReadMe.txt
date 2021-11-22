----------------------------------------

Code Descriptions:

Traveling waves in the prefrontal cortex during working memory
Sayak Bhattacharya, Scott L. Brincat, Mikael Lundqvist, Earl K. Miller
Correspondence: ekmiller@mit.edu


----------------------------------------
Main Code:


1. sb_spiral_curl_master2 - inputs LFP data and calculates circular-circular correlations

2. sb_spiral_curl_count - uses output from sb_spiral_curl_master2 and extracts the needed coefficients for further analysis - coefficient files are provided in the dropbox link below as well - this also contains a section to plot wave direction trends

3. sb_spiral_count_2D_euclidean - rotating vs planar waves from coefficients and simulated data

4. sb_spikes - analyzing spike data


-----------------------------------------

curl_planar01_sp24r_sp42r_sp24_sp42_noiseavg_fix_5curls.mat - coefficients obtained for simulated 72 types of waves . 1-32 are different directed planar waves. 33-72 are rotating waves of different directions and wavelengths (see Supplemental Figure 1).

Coefficient data files: 
Dropbox link: https://www.dropbox.com/sh/db01j12ij23ng6q/AAC2CHpQlYSJ5LLSerXplhZaa?dl=0


1. T0810_T0815_T0820_T0823_T0828_3freq_NEW_dPFC_all - coefficients for Subject 1 right hemisphere
2. T0810_T0815_T0820_T0823_T0828_3freq_NEW_dPFC_all_left - same for left hemisphere
3. T0810_T0815_T0820_T0823_T0828_3freq_NEW_dPFC_all_spikes_new - spikes for Subject 1 right hemisphere
4. T0810_T0815_T0820_T0823_T0828_3freq_NEW_dPFC_all_spikes_left_new  - same for left hemisphere

5. D0611_D0612_D0731_D0605_D0517_Dothers_3freq_NEW_dPFC_all - coefficients for Subject 2 left hemisphere
6. D0611_D0612_D0731_D0605_D0517_Dothers_3freq_NEW_dPFC_all_right - same for right hemisphere
7. D0611_D0612_D0731_D0605_D0517_Dothers_3freq_NEW_dPFC_all_spikes _new - spikes for Subject 2 left hemisphere
8. D0611_D0612_D0731_D0605_D0517_Dothers_3freq_NEW_dPFC_all_spikes_new_right - same for right hemisphere

