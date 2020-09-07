This set of scripts allow to automatically calculate Ripley's function in a set of binary masks.
Multiple sections will be aggregated by sample. 

The process consists of 2 steps:

1- RUN "1_binary_mask_ipsi.ijm" macro in Fiji (ImageJ) on each sample folder to obtain a binary mask of the ipsilateral channel

2- RUN "2_Quantify_dispersion_using_Ripley's_function.R" script on R. It will generate X-Y coordinates from the binary mask, these coordinates willd be converted to point patterns (ppp) and Ripley's function (Kest) computed.
    ![image](https://github.com/saneglaz/Image-Analysis/blob/master/Quantify_dispersion_using_Ripley's_function/Formula.jpg)

- "formula.jpg": scheme showing how distances are measured for each point; Ripley's formula
    - "explanation.jpg" explain how to interpret the results
    - "Ripley's_function.jpg" shows a cummulative function with demo results
