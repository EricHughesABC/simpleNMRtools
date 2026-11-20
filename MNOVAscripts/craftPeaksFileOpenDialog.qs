function craftPeaksFileOpenDialog() {


        // read in the json file and return the json object
    function readCRAFTpeaksFile(jsonfilename){

       var peaksTxt = "";
        var fin = new File(jsonfilename);
        if (fin.open(File.ReadOnly)) {
            var sin = new TextStream(fin, 'UTF-8');
            peaksTxt = sin.readAll();
        }
        else{
            MessageBox.warning("Could not open json file\n" + jsonfilename);
        }
        fin.close();
        return peaksTxt;
    }

    var doc = Application.mainWindow.activeDocument;
    var nmrspec = new NMRSpectrum(nmr.activeSpectrum());
    var specFileName = nmrspec.getParam("Data File Name");
    print(specFileName);

    // split specFilename to get the path upto and including craft
    var specPath = specFileName.substring(0, specFileName.lastIndexOf("craft") + 5);
    print(specPath);

    // add the rest of the path to where CRAFT output filea are located
    specPath += "/craftdbox5b/CRAFTpublish";

    var craftPeaksDialog = FileDialog.getOpenFileName("*.txt", "Read CRAFT peaks file", specPath);


    print(craftPeaksDialog);

    if (craftPeaksDialog) {
        print("Selected CRAFT peaks file: " + craftPeaksDialog);
        // read the selected CRAFT peaks file which is a txt file
        var peaksTxt = readCRAFTpeaksFile(craftPeaksDialog);
        print(peaksTxt);

        // the file consists of a header followed by peak data the data comes after the line #F1ppm     F2ppm   Amp
        // create an array of peak data lines of the format
        // peakdata[i] = { F1ppm: <value>, F2ppm: <value>, Amp: <value> }
        // convert string values to floats

        // this text can slightly vary in terms of spaces so we use a regular expression to split the columns

        // var peakDataStart = peaksTxt.indexOf("#F1ppm     F2ppm   Amp");
        var peakDataStart = peaksTxt.indexOf("#F1ppm");
        var peakDataArray = [];
        if (peakDataStart !== -1) {
            var peakData = peaksTxt.substring(peakDataStart).split("\n").slice(1);
            
            for (var i = 0; i < peakData.length; i++) {
                if (peakData[i].trim() !== "") {
                    var parts = peakData[i].trim().split(/\s+/);
                    peakDataArray.push({ F1ppm: parseFloat(parts[0]), F2ppm: parseFloat(parts[1]), Amp: parseFloat(parts[2]) });
                }
            }
        }
        print(peakDataArray);

        // create an empty mnova Peaks object
        var mnovaPeaks = nmrspec.peaks();
        for (var i = 0; i < peakDataArray.length; i++) {
            mnovaPeaks.append(Peak(peakDataArray[i].F1ppm, peakDataArray[i].F2ppm, peakDataArray[i].Amp));
        }
        nmrspec.setPeaks(mnovaPeaks);
        nmrspec.process();
    } else {
        print("No file selected.");
    }
}