input=getDirectory("Choose a Directory");//Select a folder
output=input

function action(input, output, filename) {
	
	//abrir archivo, split channels y cierra canal azul y rojo
	open(input + filename);
	run("Split Channels");
	selectWindow(filename + " (blue)");
	close();
	selectWindow(filename + " (red)");
	close();
	
	//selecciona canal ipsilateral	
	selectWindow(filename + " (green)");

	setTool("polygon");
	waitForUser("Draw a selection to calculate background and click OK");
	getStatistics(area, mean);
    run("Select None");
	slices = nSlices;
	for (i=0; i<slices; i++){
		Stack.setPosition(1, i, 1);
    	run("Subtract...", "value="+mean);
}
		
	//procesa imagen
	run("Subtract Background...", "rolling=50");
	run("Auto Threshold", "method=Moments white");
	run("RGB Color");
	
	//Save
	saveAs("PNG", output + filename);
	close();
} 

list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);

saveAs("Results", output + "Results.csv");