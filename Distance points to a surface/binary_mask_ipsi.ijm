input = getString("input?", "C:\\Users\\Frigopie\\Desktop\\DISROOT\\CUANTIFICACION SC\\EH13\\4\\");
output = getString("output?", "C:\\Users\\Frigopie\\Desktop\\DISROOT\\CUANTIFICACION SC\\EH13\\4\\");


function binary_mask(input, output, filename) {
	
	//abrir archivo, split channels y cierra canal azul y rojo
	open(input + filename);
	run("Split Channels");
	selectWindow(filename + " (blue)");
	close();
	selectWindow(filename + " (red)");
	close();
	
	//selecciona canal ipsilateral	
	selectWindow(filename + " (green)");
		
	//procesa imagen
	run("Subtract Background...", "rolling=50");
	run("Auto Threshold", "method=Moments white");
	run("RGB Color");

	saveAs("PNG", output + filename);
	close();
} 

list = getFileList(input);
for (i = 0; i < list.length; i++)
        binary_mask(input, output, list[i]);