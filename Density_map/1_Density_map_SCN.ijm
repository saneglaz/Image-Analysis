folder = getDirectory("Choose the directory with the images to analyse");

File.makeDirectory(folder + "OVERLAP");
File.makeDirectory(folder + "ROIS_FIX");
File.makeDirectory(folder + "IPSI_clean");

IPSI = folder + "IPSI_clean" + File.separator;
OIB = folder + "OIB" + File.separator;
ROIS_FIX = folder + "ROIS_FIX" + File.separator;
list = getFileList(OIB);
filestring=File.openAsString(folder + File.separator + "angles.txt");
angles=split(filestring, "\n"); 
filestring=File.openAsString(folder + File.separator + "flip.txt");
flip=split(filestring, "\n"); 
roiManager("Reset");

//START LOOP TO ANALYZE ALL IMAGES IN SELECTED FOLDER
for (i=0; i<list.length; i++){ 
	run("Set Measurements...", "area area_fraction redirect=None decimal=2");	
	//Open image and get name
	run("Bio-Formats", "open=" + OIB + list[i] + " autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	T=getTitle;
	name=substring(T, 0, lengthOf(T)-4);//Avoid .tif ending

	run("Rotate... ", "angle="+angles[i]+" grid=1 interpolation=Bilinear enlarge stack");
	if(flip[i]== "YES"){
		run("Flip Horizontally", "stack");
	}
	
	background();
	threshold_imagecalc(T);

	File.makeDirectory(IPSI + name);
	IPSI_SAMPLE =  IPSI + name + File.separator;
	File.makeDirectory(IPSI_SAMPLE + "green_side");
	File.makeDirectory(IPSI_SAMPLE + "red_side");
	GREEN_SIDE = IPSI_SAMPLE + "green_side" + File.separator;
	RED_SIDE = IPSI_SAMPLE + "red_side" + File.separator;
	
	selectWindow("red_");
	run("Duplicate...", "duplicate");
	rename("647_");
	selectWindow("green_");
	run("Duplicate...", "duplicate");
	rename("488_");
	
	open(ROIS_FIX + name + "_green.roi"); 
	if(getTitle()== name + "_green.roi"){
		roiManager("Add");
		close();
		roiManager("select", 0);
	}
	roiManager("Add");
//	waitForUser("greenside verde?");
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "stack");	
	run("Make Inverse");
	run("Crop");
	run("Canvas Size...", "width=679 height=830 position=Center zero");
	run("Image Sequence... ", "format=PNG save=" + GREEN_SIDE + T + "_0000.png");
	close();
	
	selectWindow("647_");
	roiManager("select", 0);
//	waitForUser("greenside rojo?");	
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "stack");	
	run("Make Inverse");
	run("Crop");
	run("Canvas Size...", "width=679 height=830 position=Center zero");
	run("Image Sequence... ", "format=PNG save=" + GREEN_SIDE + T + "_0000.png");
	close();

	selectWindow("red_");
	rename("647_");
	selectWindow("green_");
	rename("488_");
	
	open(ROIS_FIX + name + "_red.roi"); 
	if(getTitle()== name + "_red.roi"){
		roiManager("Add");
		close();
		roiManager("select", 1);
	}
	roiManager("Add");
//	waitForUser("redside verde?");
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "stack");	
	run("Make Inverse");
	run("Crop");
	run("Flip Horizontally", "stack");
	run("Canvas Size...", "width=679 height=830 position=Center zero");
	run("Image Sequence... ", "format=PNG save=" + RED_SIDE + T + "_0000.png");
	close();
	
	selectWindow("647_");
	roiManager("select", 1);
//	waitForUser("redside rojo?");	
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "stack");	
	run("Make Inverse");
	run("Crop");
	run("Flip Horizontally", "stack");
	run("Canvas Size...", "width=679 height=830 position=Center zero");
	run("Image Sequence... ", "format=PNG save=" + RED_SIDE + T + "_0000.png");

	//close();
	roiManager("Reset");
	run("Close All");	
}

//------------------------------
//FUNCTIONS

function threshold_imagecalc(T) {
	run("Split Channels");
	selectWindow("C1-"+T);
	run("Enhance Contrast...", "saturated=0.1 normalize process_all");
	run("Convert to Mask", "method=IsoData background=Default black");
	invert_bug();

	selectWindow("C3-"+T);
	run("Enhance Contrast...", "saturated=0.1 normalize process_all");
	run("Convert to Mask", "method=IsoData background=Default black");
	invert_bug();

	selectWindow("C2-"+T);
	run("Enhance Contrast...", "saturated=0.1 normalize process_all");
	run("Convert to Mask", "method=Moments background=Dark calculate black");
	//run("Convert to Mask", "method=Momments black");
	rename("af10");
	//invert_bug();

	waitForUser("COMPROBAR IMAGENES ANTES COLOC");

	imageCalculator("AND create stack", "C1-"+T,"af10");
	imageCalculator("AND create stack", "C3-"+T,"af10");
	

 	imageCalculator("Subtract create stack", "Result of " + "C1-"+T,"Result of " + "C3-"+T);
	rename("green_");
	invert_bug();

	imageCalculator("Subtract create stack", "Result of " + "C3-"+T,"Result of " + "C1-"+T);
	rename("red_");
	invert_bug();

	waitForUser("COMPROBAR IMAGENES DESPUES COLOC");
	selectWindow("C1-"+T);
	close();
	selectWindow("C3-"+T);
	close();
	selectWindow("af10");
	close();
}

function background() {

	Stack.setPosition(2,15, 1);
	setTool("polygon");
	waitForUser("Draw a selection to calculate background and click OK");
	getStatistics(area, mean);
    run("Select None");
	slices = nSlices;
	for (i=0; i<slices; i++){
		Stack.setPosition(2, i, 1);
    	run("Subtract...", "value="+mean);
}
	Stack.setPosition(2,15, 1);
	run("Subtract Background...", "rolling=10 stack");
	Stack.setPosition(1,15, 1);
	run("Subtract Background...", "rolling=10 stack");
	Stack.setPosition(3,15, 1);
	run("Subtract Background...", "rolling=10 stack");	

}

function invert_bug() {
	run("Clear Results");
    run("Measure");
	area_fraction = getResultString("%Area");
	if(area_fraction> 60){
		run("Invert");
	}}