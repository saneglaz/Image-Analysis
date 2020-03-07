This set of scripts allow to automatically calculate in a set of images the distances of positive pixels to an user defined segmented line.
Multiple sections will be aggregated by sample. 
The process consists of 4 steps:

1- RUN "1_binary_mask_ipsi.ijm" macro in Fiji (ImageJ) on each sample folder to obtain a binary mask of the secondary channel
   - The user must draw the line in first place

2- RUN "2_4_Distance_points_to_a_surface.R" script on R only up to section "#Generate files with the coordinates" included. It will generate X-Y coordinates from the binary mask.

3- RUN "3_Distance_points_to_a_surface.ijm"  macro in Fiji (ImageJ) on each sample folder to obtain a binary mask of the secondary channel. It will read X-Y coordinates to create points and calculate their distance to the input line.

4- RUN "2_4_Distance_points_to_a_surface.R" from "#Read distances" section. 
   - It will read the distance data to generate vectors for each sample by aggregating sections. 
   - Vectors are later reescaled in order to have the same number of points per sample. 
   - Statistic analyses include:
      * Cummulative distribution for each group. The envelope of each function represents the SEM (standar error of the mean).
      * Violin plot

