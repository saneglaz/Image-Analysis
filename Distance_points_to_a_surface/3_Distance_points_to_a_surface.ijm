input=getDirectory("Choose a Directory");//Select a folder
output=input

PNG = input + "PNG\\";
PUNTOS = input + "PUNTOS\\";
ROIS = input + "ROIS\\";

function rois(input, ROIS, filename) {
	open(PNG + filename); 
 	name = getTitle();
 	roiManager("Show All");
 	roiManager("Reset");
 	run("Clear Results");
 	run("Flip Horizontally");
 	run("Rotate 90 Degrees Left");
        
    // Draw Curves

	setTool("polyline");
	waitForUser("Draw your freehand curve and click OK");
	run("Fit Spline");
	saveAs("Selection", ROIS + name + ".png");

	close();       
} 	


function distance(input, PNG, filename) {
	open(input + "PNG\\" + filename);  
 	name = filename;
 	roiManager("Show All");
 	roiManager("Reset");
 	run("Clear Results");
 	run("Flip Horizontally");
 	run("Rotate 90 Degrees Left");
	
	open(input + "ROIS\\" + filename + ".roi");
 	Roi.setName("Curve");
	roiManager("Add");

 	fileName = input + "PUNTOS\\" + filename + "_points.csv";
	allText = File.openAsString(fileName);
	tmp = split(fileName,".");
	// get file format {txt, csv}
	posix = tmp[lengthOf(tmp)-1];
	// parse text by lines
	text = split(allText, "\n");
 
	// define array for points
	var xpoints = newArray;
	var ypoints = newArray; 
 
	// in case input is in TXT format
	if (posix=="txt") {	
		print("importing TXT point set...");
		//these are the column indexes
		hdr = split(text[0]);
		nbPoints = split(text[1]);
		iX = 0; iY = 1;
		// loading and parsing each line
		for (i = 2; i < (text.length); i++){
	   		line = split(text[i]," ");
	   		setOption("ExpandableArrays", true);   
	   		xpoints[i-2] = parseInt(line[iX]);
	   		ypoints[i-2] = parseInt(line[iY]); 
	   		print("p("+i-1+") ["+xpoints[i-2]+"; "+ypoints[i-2]+"]"); 
		} 	
		// in case input is in CSV format
	} else if (posix=="csv") {
		print("importing CSV point set...");
		//these are the column indexes
		hdr = split(text[0]);
		iLabel = 0; iX = 1; iY = 2;
		// loading and parsing each line
		for (i = 1; i < (text.length); i++){
	   		line = split(text[i],",");
	   		setOption("ExpandableArrays", true);   
	   		xpoints[i-1] = parseInt(line[iX]);
	   		ypoints[i-1] = parseInt(line[iY]);
	   		print("p("+i+") ["+xpoints[i-1]+"; "+ypoints[i-1]+"]"); 
			} 
	// in case of any other format
	} else {
		print("not supported format...");	
	}
 
	// show the points in the image
	makeSelection("point", xpoints, ypoints); 

	Roi.setName("Points");
	roiManager("Add");

	// 3. Make Distance Map
	getDimensions(x,y,c,z,t);
	newImage(name+" - Distance Map", "8-bit white", x, y, 1);
	setForegroundColor(0, 0, 0);
	roiManager("Select",0);
	run("Draw", "slice");
	setAutoThreshold("Default dark");
	run("Convert to Mask");
	run("Exact Euclidean Distance Transform (3D)");

	// 4. Measure Points
	roiManager("Select", 1);
	run("Measure");

// Result appears in "Mean Column
	
	saveAs("Results", output + filename + "_distances.csv");

	
	selectWindow(filename + " - Distance Map");
	close();
	selectWindow("EDT");
	close();
	selectWindow(filename);
	close();

}


list=getFileList(PNG); 
for (i = 0; i < list.length; i++) {
   rois(input, ROIS, list[i]);
} 

list=getFileList(PNG); 
for (i = 0; i < list.length; i++) {
   distance(input, ROIS, list[i]);
} 
