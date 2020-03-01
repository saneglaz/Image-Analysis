RUN "R_distribution.ijm" macro on Fiji (ImageJ) for each sample folder

- The macro will process the images for each section. Input TIFF/PNG files should be placed in "PNG" folder. 
- Roi selections can be previously placed in "ROI" folder
- The loop will check if the roi selection for the file exist. If not, it will be required to the user.
- OUTPUT: In the "R_DISTRIBUTION" folder, one file for each flurescence channel will be created ("filename_contra.txt" & "filename_ipsi.txt") 
-FUNCTIONS:
	- background: Allow to select an empty region to calculate the backgroud, and remove that intensity value from the whole image
	- textimg: Preprocess each channel and obtain an intensity matrix for each channel
	- LUT: Generates a pseudocolor image for the R distribution in JPEG format. It should be used after generate the R-distribution files with the below R-script.


USE "R_distribution.R" script in R to calculate the R distribution and analyze the data:

- The loop will navigate through the folder tree processing all samples
- OUTPUT: 	-"filename_Rdist.txt" file contains a matrix with the R-values in the sample
		-results data.frame
- The resulting data will be integrated in the "results" data.frame including:
	- ratio_overlap: ratio of overlaping area: (R[-0.5, 0.5] / total area)
	- ratio_ipsidom: ratio of the pixels enriched for the ipsilateral channel (R[1.75, 2.5] / total area)
	- ratio_contradom: ratio of the pixels enriched for the contralateral channel (R [-2.5, -1.75] / total area)
	- total_px: total area
- To better visualize the results, the code to generate an histogram is included: The vectors with the R-values from each sample will be reescaled and all the samples for each group averaged  

-FUNCTIONS:
	- R_distribution: open the "filename_contra.txt" & "filename_ipsi.txt" and compute the R distribution
		* The intensity threshold was set at 5. Lower intensities were considered as "empty".
		* "Empty" pixels in both channels were discarded
		* R=log10(Fi/Fc); Fi is the fluorescence intensity in the ipsilateral channel and Fc is the fluorescence intensity in the contralateral channelç
		* To allow the calculation of the logarithm, the remaining zero-value pixels were replaced by 0.01
	- se: calculate the standard error of the mean

