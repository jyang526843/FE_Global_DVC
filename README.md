# Finite element Global Digital Volume Correlation (FE_Global_DVC)
Finite-element-based global DVC method (guarantee global kinematic compatibility and decrease noise by adding regularization penalties).  
 
## Prerequisites & Installation
FE_Global_DVC MATLAB code was tested on MATLAB versions later than R2018a. Please download and unzip the code to the MATLAB working path. Then, execute the mail file main_FE_Global_DIC.m.

## Code Manual
It is already included in the code. It's also available at my Researchgate: https://www.researchgate.net/publication/xxx

## FE-Global-DVC example dataset
Example images download link: 
https://uwmadison.box.com/s/kr6pjje9yfi9k9va6cbk1ayvhvbs5nvq <p>
  
## ****** ATTENTION ******  
% The "x,y,z" or "1-,2-,3-" coordinates in the ALDVC code always correspond to the 1st, 2nd and 3rd indices of Matlab workspace variable. For example, p_meas(:,1) and p_meas(:,2) are the x- & y-coordinates of scattered points.  
 
% This is a little different from some MATLAB image processing functions. 
% For example, if a 3D image has size MxNxL, in this code, we always have the image size_x=M, size_y=N, size_z=L. If you use some Matlab computer vision/image post-processing function, for example, 'imagesc3D', or 'imshow3D', or 'surf', it will reads size_x=N, size_y=M, size_z=L. 
 
% Please pay attention to this difference.  

## Citation
If used please cite
```bibtex
@article{Yang2020aldvc,
  title={Augmented Lagrangian Digital Volume Correlation (ALDVC)},
  author={Yang, J. and Hazlett, L. and Landauer, A. K. and Franck, C.},
  journal={Experimental Mechanics},
  year={2020},
  Url={https://doi.org/10.1007/s11340-020-00607-3}
}
```
 
[1] J Yang, L Hazlett, AK Landauer, C Franck. Augmented Lagrangian Digital Volume Correlation (ALDVC). Experimental Mechanics, 2020.  

Full text can be requested at: 
* Exp Mech Website: https://link.springer.com/article/10.1007/s11340-020-00607-3
* ResearchGate: https://www.researchgate.net/publication/343676441_Augmented_Lagrangian_Digital_Volume_Correlation_ALDVC

## Contact and support
Email: aldicdvc@gmail.com;  -or- Jin Yang, jyang526@wisc.edu; -or- Prof. Christian Franck, cfranck@wisc.edu


##
 
<p align="center">
  <img width="538" height="301" src="https://github.com/FranckLab/ALDVC/blob/master/aldvc_logo.png">
</p>

