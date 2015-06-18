# LV-segmentation-in-cardiac-MRI
semi-automatic segmentation of LV in cardiac MR images

This code performs semi-automatic segmentation of heart chambers, mainly for LV, in cardiac MR images.
Semi-automatic means that the initial contour should be given to the method by the user.

The method employes an active contour method. A shape can be incorporated into the method to resolve some of the 
drawbacks of active contour methods such as leakage and shrinkage. Defual value for the weight of the shape is zero.

Some sample MR images are provided here. 
If your the code, an image will be shown to the user and asked for an initial contour.
The initial contour can be a simple circle or any sort of delineations around or inside the LV chamber.

After initialization by the user, the automatic segmentation will be performed.
The final figure will show a comparison between the automatic and ground truth segmentations. The ground truth is provided
by expert radiologists.  

References 
1- Chan, T.F.; Vese, L.A., "Active contours without edges," Image Processing, IEEE Transactions on , vol.10, no.2, pp.266,277, Feb 2001
doi: 10.1109/83.902291

2- Pluempitiwiriyawej, C.; Moura, J.M.F.; Yi-Jen Lin Wu; Chien Ho, "STACS: new active contour scheme for cardiac MR image segmentation," Medical Imaging, IEEE Transactions on , vol.24, no.5, pp.593,603, May 2005
doi: 10.1109/TMI.2005.843740

