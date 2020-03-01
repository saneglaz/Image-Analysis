folder = getDirectory("Choose the directory with the images to analyse");

File.makeDirectory(folder + "R_DISTRIBUTION");
PNG = folder + "PNG" + File.separator;
ROIS = folder + "ROIS" + File.separator;
R_DISTRIBUTION = folder + "R_DISTRIBUTION" + File.separator;

run("Clear Results"); 
list = getFileList(PNG);
arrayNames = newArray(list.length);

//START LOOP TO ANALYZE ALL IMAGES IN SELECTED FOLDER
for (i=0; i<list.length; i++){ 
	//Open and get data
	open(PNG + list[i]);
	filename=getTitle();
	name=substring(filename, 0, lengthOf(filename)-4);//Avoid .tif ending
	arrayNames[i] = name;

	if(File.exists(R_DISTRIBUTION + name + "_contra.txt")==0){
		textimg(R_DISTRIBUTION,filename,name);
		waitForUser("Process data in R", "Thereafter press OK");	
		if(File.exists(R_DISTRIBUTION + name + "_Rdist.txt")==1){
		LUT(R_DISTRIBUTION,name);}}
	else{
		if(File.exists(R_DISTRIBUTION + name + "_Rdist.txt")==1){
		LUT(R_DISTRIBUTION,name);}}
	run("Close All");
	}
	roiManager("Deselect");
	roiManager("Delete");
	
//------------------------------
//FUNCTIONS
	
function textimg(R_DISTRIBUTION, filename, name) {
	run("Split Channels");
	selectWindow(filename + " (blue)");
	close();

	selectWindow(filename + " (green)");
	background();
	run("Enhance Contrast...", "saturated=0.3 normalize");
	
	 
	open(ROIS + filename + "_contra.roi"); 
	if(getTitle()==filename+"_contra.roi"){
		roiManager("Add");
		close();
		roiManager("select", 0);
	}
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "slice");	
	
	saveAs("Text Image", R_DISTRIBUTION + name + "_contra.txt");
	selectWindow(filename + " (green)");
	close();

	selectWindow(filename + " (red)");
	background();
	run("Enhance Contrast...", "saturated=0.3 normalize");
	
	open(ROIS + filename + "_contra.roi"); 
	if(getTitle()==filename+"_contra.roi"){
		roiManager("Add");
		selectWindow(filename + "_contra.roi");
		close();
		roiManager("select", 0);
	}
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "slice");
	saveAs("Text Image", R_DISTRIBUTION + name  + "_ipsi.txt");

	selectWindow(filename + " (red)");
	close();
}

function background() {
	setTool("polygon");
	waitForUser("Draw a selection to calculate background and click OK");
          getStatistics(area, mean);
          run("Select None");
          run("Subtract...", "value="+mean);
      }

function LUT(R_DISTRIBUTION, name) {

	run("Text Image... ", "open="+R_DISTRIBUTION + name+ "_Rdist.txt");
	run("32_Colors");
	run("Invert LUT");


	saveAs("Jpeg", R_DISTRIBUTION + name + ".Rdist.jpg");
	close();
}

        