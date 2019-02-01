# SPIN-VM
Spatially Informed Voxelwise Modeling

This technique is described in the [following paper](https://www.sciencedirect.com/science/article/pii/S1053811918321256):
* Çelik, E., Dar, S. U. H., Yılmaz, Ö., Keleş, Ü., & Çukur, T. (2019). Spatially informed voxelwise modeling for naturalistic fMRI experiments. NeuroImage, 186, 741-757. doi:10.1016/j.neuroimage.2018.11.044

Relevant files can be found under the **data** folder. semS1_presp_R.mat stores the preprocessed BOLD response file for a single ROI (PPA). Response file should be of size (time points X number of voxels) and stimulus file should be of size (time points X number of features). semgab2_movie_stim.mat stores the preprocessed stimulus file. ROI_S1.mat stores ROI indices, and semS1_Rvallcvall_noise_cf.hf5 stores noise correction factors. Run spinvm_demo.m for a quick demo.

You are encouraged to modify/distribute this code. However, please acknowledge this code and cite the paper appropriately. 

For any questions, comments and contributions, please contact Tolga Cukur (cukur[at]ee.bilkent.edu.tr)

(c) ICON Lab 2019
