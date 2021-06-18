// Quant_Biofilms.ijm - ImageJ macro to quantify biofilm images in multi-series file
//
// Usage: macro prompts user to select multi-series file containing images to analyze
// Results: new Results table with 1 row of measurements per input image
// Requirements: Bio-Formats plugin to open multi-series file (recommend Fiji)
//
// Copyright: Dundee Imaging Facility, University of Dundee (2021)
// Author: Graeme Ball (g.ball@dundee.ac.uk)
// License: GNU GPL
//
// Version 6
// v1: iterate over multi-series file, find colony, measure area & intensities
//     show results, save output folder of snapshot images
// v2: add Pearsons and Manders M1, M2; save results csv to output folder
// v3: dialog to adjust parameters & make Manders coefficients optional
//     - support single fluorescent channel (no coloc stats!)
//     - do not round percentages; update results during macro
// v4: add optional "subcolony" segmentation and test options
//     - "segment subcolonies" option where >1 colony of interest
//     - option to specify firstSeries and lastSeries to test on large datasets
//     - hard-coded testSegmentation parameter (=false) and function
// v5: add greenForegroundTotal and redForegroundTotal stats & accept BF-only
//     - total intensity of foreground pixels (N.B. background *not* subtracted)
//     - colonyArea for brightfield-only images (green AND red channelNumber=0)
//     - fix bug: subcolony area was pixels, now calibrated units
//     - optional background subtraction for non-clear plates
// v6: add test segmentation option to user dialog
//     - fix bug: month was 0-based in output folder time stamp, now 1-based
//

// -- parameters & options --
channelNumberBF = 1;                 // brightfield
channelNumberGreen = 2;
channelNumberRed = 3;
minimumColonyArea = 5000000;         // minimum area (in um^2)
backgroundGreen = 250;               // 113 = colony mean+5sd from green -ve, 5852
backgroundRed = 150;                 // 18 = colony mean+5sd from red -ve, 1473
displayMaxGreen = 4095;
displayMaxRed = 4095;
doSBg = false;                        // do "Subtract Background" (intended for 1-channel milk)
rbRad = 100;                          // rolling ball radius for Subtract Background
segmentSubcolonies = false;           // segmentation method includes small subcolonies
doCalculateManders = false;
specifyFirstAndLast = false;         // specify first & last image to analyze
testSegmentation = false;            // display segmentation results only
firstSeries = 0;
lastSeries = 1;
setOption("BlackBackground", true);  // white objects, black background


// -- Main Macro --

// Clean up before starting
run("Close All");
if (isOpen("ROI Manager")) {
    selectWindow("ROI Manager");
    run("Close");
}
run("Clear Results");

// user update of parameters
Dialog.create("Quant_Biofilms");
Dialog.addNumber("Brightfield channel number", channelNumberBF);
Dialog.addNumber("Green channel number (0=none)", channelNumberGreen);
Dialog.addNumber("Red channel number (0=none)", channelNumberRed);
Dialog.addNumber("Minimum (sub)colony area (unit^2)", minimumColonyArea);
Dialog.addNumber("Green channel background level", backgroundGreen);
Dialog.addNumber("Red channel background level", backgroundRed);
Dialog.addNumber("Green display max", displayMaxGreen);
Dialog.addNumber("Red display max", displayMaxRed);
Dialog.addCheckbox("Subtract Background?", doSBg);
Dialog.addCheckbox("Segment subcolonies?", segmentSubcolonies);
Dialog.addCheckbox("Calculate Manders coefficients?", doCalculateManders);
Dialog.addCheckbox("Specify first & last image to analyze?", specifyFirstAndLast);
Dialog.addCheckbox("Test segmentation only?", testSegmentation);

Dialog.show();
channelNumberBF = Dialog.getNumber();
channelNumberGreen = Dialog.getNumber();
channelNumberRed = Dialog.getNumber();
minimumColonyArea = Dialog.getNumber();
backgroundGreen = Dialog.getNumber();
backgroundRed = Dialog.getNumber();
displayMaxGreen = Dialog.getNumber();
displayMaxRed = Dialog.getNumber();
doSBg = Dialog.getCheckbox();
segmentSubcolonies = Dialog.getCheckbox();
doCalculateManders = Dialog.getCheckbox();
specifyFirstAndLast = Dialog.getCheckbox();
testSegmentation = Dialog.getCheckbox();

if (specifyFirstAndLast) {
	Dialog.create("Quant_Biofilms");
	Dialog.addNumber("First image to analyze", firstSeries);
	Dialog.addNumber("Last image to analyze", lastSeries);
	Dialog.show();
	firstSeries = Dialog.getNumber();
	lastSeries = Dialog.getNumber();
}

// determine which channels to show in merge snapshot
if (channelNumberGreen > 0 && channelNumberRed > 0) {
	// show both color channels but not brightfield
	active = newArray(1, 1, 1);
	active[channelNumberBF-1] = 0;
	snapshotActiveChannels = "" + active[0] + "" + active[1] + "" + active[2];
} else if (channelNumberGreen > 0 || channelNumberRed > 0) {
	// assume brightfield & color, show both
	snapshotActiveChannels = "11";
} else {
	// 1 channel only (brightfield)
	snapshotActiveChannels = "1";
}
numberOfChannelsExpected = snapshotActiveChannels.length;

// Open multi-series file
msFile = File.openDialog("Choose multi-series file containing biofilm images to analyze");
inputFolder = File.getParent(msFile);
msFilename = File.getName(msFile);
if (!testSegmentation) {
	outputFolder = inputFolder + File.separator + baseName(msFilename) + "_quantBF_" + timeStamp();
	File.makeDirectory(outputFolder);
}
run("Bio-Formats Macro Extensions");
Ext.setId(msFile);  // initialize file
Ext.getSeriesCount(seriesCount);
// set firstSeries and lastSeries if not specified
if (!specifyFirstAndLast) {
	firstSeries = 0;
	lastSeries = seriesCount - 1;
}
setBatchMode("hide");
run("Set Measurements...", "  redirect=None decimal=3");  // to avoid AnalyzeParticles junk
// Iterate over image series in file writing 1 result line per series
for (s = firstSeries; s <= lastSeries; s++) {
	showProgress((s-firstSeries) / (lastSeries-firstSeries+1));
	Ext.setSeries(s);
    Ext.getSeriesName(imageName);
    Ext.getSizeC(nChannels);
    if (nChannels != numberOfChannelsExpected) {
    	message = "skipping image " + imageName + " - expected ";
    	message = message + numberOfChannelsExpected + " channels, found " + nChannels; 
    	print(message);
    } else if (testSegmentation) {
    	Ext.openImagePlus(msFile);
		imageID = getImageID();
		if (segmentSubcolonies) {
			colonyROIindex = addSubcoloniesRoiToManager(channelNumberBF, minimumColonyArea, doSBg, rbRad);
		} else {
			colonyROIindex = addColonyRoiToManager(channelNumberBF, minimumColonyArea, doSBg, rbRad);
		}
		selectImage(imageID);
		if (colonyROIindex > -1) {
			roiManager("select", colonyROIindex);
			run("Add Selection...");  // add to overlay
		}
		run("Enhance Contrast", "saturated=0.35");
		roiManager("reset");
    } else {
	    Ext.openImagePlus(msFile);
		Stack.getUnits(xu, yu, zu, tu, vu);
		areaUnits = "_" + xu + "2";  // probably micron^2
		updateResults();  // otherwise 1st "Label" entry goes missing...
		row = nResults;
		imageName = getTitle();  // sometimes empty from Ext.getSeriesName?
		setResult("Label", row, imageName);
		// Find colony, measure & report stats
		if (segmentSubcolonies) {
			colonyROIindex = addSubcoloniesRoiToManager(channelNumberBF, minimumColonyArea, doSBg, rbRad);
		} else {
			colonyROIindex = addColonyRoiToManager(channelNumberBF, minimumColonyArea, doSBg, rbRad);
		}
		colonyArea = 0;
		greenMean = 0;
		greenMax = 0;
		greenStd = 0;
		greenTotal = 0;
		redMean = 0;
		redMax = 0;
		redStd = 0;
		redTotal = 0;
		if (colonyROIindex > -1) {
			// colony stats
			roiManager("select", colonyROIindex);
			roiManager("rename", "colony_" + safeName(imageName));
			run("Add Selection...");  // add to overlay
			getStatistics(colonyArea);
			setResult("colonyArea" + areaUnits, row, colonyArea);
			// green channel stats
			if (channelNumberGreen > 0) {
				Stack.setChannel(channelNumberGreen);
				getRawStatistics(nPixels, greenMean, minG, greenMax, greenStd);
				greenTotal = nPixels * greenMean;
				I_G = getPixelIntensities();
				pctGrn = calcPercentAboveThresh(backgroundGreen);
				greenForegroundTotal = calcForegroundTotal(I_G, backgroundGreen);
				setResult("greenMean", row, greenMean);
				setResult("greenMax", row, greenMax);
				setResult("greenStd", row, greenStd);
				setResult("greenTotal", row, greenTotal);
				setResult("greenForegroundTotal", row, greenForegroundTotal);
				setResult("pctGreen", row, pctGrn);
			}
			// red channel stats
			if (channelNumberRed > 0) {
				Stack.setChannel(channelNumberRed);
				getRawStatistics(nPixels, redMean, minR, redMax, redStd);
				redTotal = nPixels * redMean;
				I_R = getPixelIntensities();
				pctRed = calcPercentAboveThresh(backgroundRed);
				redForegroundTotal = calcForegroundTotal(I_R, backgroundRed);
				setResult("redMean", row, redMean);
				setResult("redMax", row, redMax);
				setResult("redStd", row, redStd);
				setResult("redTotal", row, redTotal);
				setResult("redForegroundTotal", row, redForegroundTotal);
				setResult("pctRed", row, pctRed);
			}
			// coloc stats
			if (channelNumberGreen > 0 && channelNumberRed > 0) {
				PCC = calcPCCrho(I_G, I_R);
				threshPixCountAND = thresholdIntensities2C_AND(I_G, I_R, backgroundGreen, backgroundRed);
				I_Gand = Array.trim(I_G, threshPixCountAND);
				I_Rand = Array.trim(I_R, threshPixCountAND);
				PCC_obj = calcPCCrho(I_Gand, I_Rand);
				setResult("PCC", row, PCC);
				setResult("PCC_obj", row, PCC_obj);
				if (doCalculateManders) {
					M1M2 = calcManders(I_R, backgroundRed, I_G, backgroundGreen);
					M1 = M1M2[0];
					M2 = M1M2[1];
					setResult("M1_R", row, M1);
					setResult("M2_G", row, M2);
				}
			}
		}
		// Save snaphot image with fixed display settings & colony ROI
		if (nChannels > 1) {
			Stack.setChannel(channelNumberBF);
		}
		run("Grays");
		run("Enhance Contrast", "saturated=0.35");
		if (channelNumberGreen > 0) {
			Stack.setChannel(channelNumberGreen);
			run("Green");
			setMinAndMax(backgroundGreen, displayMaxGreen);
		}
		if (channelNumberRed > 0) {
			Stack.setChannel(channelNumberRed);
			run("Red");
			setMinAndMax(backgroundRed, displayMaxRed);
		}
		if (nChannels > 1) {
			Stack.setDisplayMode("composite");
			Stack.setActiveChannels(snapshotActiveChannels);
			run("RGB Color");  // should copy overlay
		}
		saveAs(outputFolder + File.separator + safeName(imageName) + "_QBf.tif");
		close();  // RGB/single-channel
		if (nChannels > 1) {
			close();  // composite
		}
		updateResults();
		roiManager("reset");
    }
}
showProgress(1);
setBatchMode("exit and display");
if (testSegmentation) {
	showMessage("Test Segmentation Results Shown");
} else {
	Table.rename("Results", "Results_Quant_Biofilms");
	saveAs("Results", outputFolder + File.separator + "Results_Quant_Biofilms.csv");
	run("Set Measurements...", "area mean standard min redirect=None decimal=3");  // reset to basics
}


// -- function definitions --

function baseName(filename) {
    // return filename string without extension
    return substring(filename, 0, lastIndexOf(filename, "."));
}

function safeName(seriesName) {
	// return series name string with filesystem safe characters
	result = replace(seriesName, "[]\\[\\(\\)/\\\\]+", "_");
	result = replace(result, "[@$%+#^&]+", "-");
	result = replace(result, "[ ]+", "_");
	return result;
}

function timeStamp(){
    // generate a time stamp string
    getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
    timeString = toString(year) + "-" + twoDigit(month+1) + "-" + twoDigit(dayOfMonth);
    DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
    timeString = timeString + "_" + DayNames[dayOfWeek];
    timeString = timeString + twoDigit(hour) + "-" + twoDigit(minute) + "-" + twoDigit(second);
    return timeString;
}

function twoDigit(n) {
    // return number as 2-digit string
    return IJ.pad(n, 2);
}

function addColonyRoiToManager(channelNumber, minimumColonyArea, doSBg, rbRad) {
	// find colony in image by autothreshold, fill, wand at center
	// add ROI to manager and return ROI index, or -1 if not found
	autoThresholdMethod = "Triangle dark";  // 	
	colonyROIindex = -1;
	run("Select None");
	getDimensions(w, h, nc, nz, nt);
	if (nc > 1) {
		Stack.setChannel(channelNumberBF);
	}
	run("Duplicate...", " ");
	if (doSBg) {
		run("Subtract Background...", "rolling=" + rbRad + " sliding");
	}
	setAutoThreshold(autoThresholdMethod);
	run("Convert to Mask");
	run("Fill Holes");
	doWand(getWidth()/2, getHeight()/2, 0, "Legacy");
	getStatistics(area);
	if (selectionType() > 0 && area > minimumColonyArea) {
		roiManager("Add");
		colonyROIindex = roiManager("count")-1;
	}		
	close();  // close mask image
	return colonyROIindex;
}

function addSubcoloniesRoiToManager(channelNumber, minimumColonyArea, doSBg, rbRad) {
	// find all subcolonies in image by autothreshold, fill, minimum size
	// add ROI to manager and return ROI index, or -1 if not found
	autoThresholdMethod = "Triangle dark";  // 	
	colonyROIindex = -1;
	inputID = getImageID();
	run("Select None");
	getDimensions(w, h, nc, nz, nt);
	if (nc > 1) {
		Stack.setChannel(channelNumberBF);
	}
	run("Duplicate...", " ");
	bfID = getImageID();
	if (doSBg) {
		run("Subtract Background...", "rolling=" + rbRad + " sliding");
	}
	setAutoThreshold(autoThresholdMethod);
	run("Create Mask");
	maskID = getImageID();
	run("Analyze Particles...", "size=" + minimumColonyArea + "-Infinity circularity=0.1-1.00 add slice");
	if (roiManager("count") > 1) {
		roiManager("deselect");
		roiManager("Combine");
		roiManager("add");
		colonyROIindex = roiManager("count")-1;
	} else if (roiManager("count") > 0) {
		colonyROIindex = roiManager("count") - 1;
	}
	selectImage(bfID);
	close();
	selectImage(maskID);
	close();
	selectImage(inputID);
	return colonyROIindex;
}

function getPixelIntensities() {
	// return an array of pixel intensities for the current image slice
	// - only within active selection if there is one

	if (selectionType() > -1) {
		// return pixel intensities only within active selection
		Roi.getContainedPoints(xc, yc);
		I = newArray(xc.length);
		for (i = 0; i < xc.length; i++) {
			I[i] = getPixel(xc[i], yc[i]);
		}
	} else {
		w = getWidth;
		h = getHeight;
		I = newArray(w * h);
		i = 0;
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				I[i] = getPixel(x, y);
				i++;
			}
		}
	}
	return I;
}

function calcPercentAboveThresh(thresh) {
	// percentage pixels above threshold for active channel & selection
	hadSelection = false;
	if (selectionType() > -1) {
		// return pixel intensities only within active selection
		hadSelection = true;
		Roi.getContainedPoints(xc, yc);
	} else {
		getDimensions(w, h, nc, nz, nt);
		makeRectangle(0, 0, w, h);
		Roi.getContainedPoints(xc, yc);
	}
	nPixels = xc.length;
	nThresholded = 0;
	for (i = 0; i < nPixels; i++) {
		 if (getPixel(xc[i], yc[i]) >= thresh) {
		 	nThresholded++;
		 }
	}
	if (!hadSelection) {
		run("Select None");
	}
	return 100 * nThresholded / nPixels;
}

function calcForegroundTotal(I, background) {
	// return total foreground intensity for pixels in I
	foregroundTotal = 0;
	for (i = 0; i < I.length; i++) {
		if (I[i] > background) {
			foregroundTotal += I[i];
		}
	}
	return foregroundTotal;
}

function thresholdIntensities2C_AND(I1, I2, thresh1, thresh2) {
	// modify arrays I1, I2 in-place keeping only above-threshold intensities
	// checks whether pixel intensities >= thresh in BOTH channels
	// return count of valid elements in I1, I2 after threhsolding (need to trim!)
	count = 0;
	for (i = 0; i < I1.length; i++) {
		if ((I1[i] >= thresh1) && (I2[i] >= thresh2)) {
			I1[count] = I1[i];
			I2[count] = I2[i];
			count++;
		}
	}
	I1 = Array.trim(I1, count);
	I2 = Array.trim(I2, count);
	return count;
}

function calcPCCrho(i1, i2) {
	// return Pearson's correlation coefficient rho
	Array.getStatistics(i1, min1, max1, mean1, std1);
	Array.getStatistics(i2, min2, max2, mean2, std2);
	rho = 0;
	for (i = 0; i < i1.length; i++) {
		rho += (i1[i] - mean1) * (i2[i] - mean2);
	}
	rho = rho / (std1 * std2 * (i1.length - 1));
	return rho;
}

function calcManders(I1, threshR, I2, threshG) {
	// calculate and report thresholded Manders coeficients M1 and M2
	// where for each pixel I1=red intensities, I2=green intensities
	// return array [M1, M2]
	//   M1 = Sum_i(Ri,coloc) / Sum_i(Ri)
	//   M2 = Sum_i(Gi,coloc) / Sum_i(Gi)
	// where:
	//   - only above-threshold intensities are considered
	//     i.e. pre-calc subtraction of threshR and threshG
	//   - Ri colocalized if Gi > O and Gi colocalized if Ri > 0
	// see: Manders et al. 1992, J. Micros.
	R = Array.copy(I1);  // channel 1, "red" intensities
	for (i = 0; i < R.length; i++) {
		R[i] = maxOf(R[i] - threshR, 0);  // subtract threshR, minimum result=0
	}
	G = Array.copy(I2);  // channel 2, "green" intensities
	for (i = 0; i < G.length; i++) {
		G[i] = maxOf(G[i] - threshG, 0);  // subtract threshG, minimum result=0
	}
	// calculate M1
	sum_Ri_coloc = 0;
	sum_Ri = 0;
	for (i = 0; i < R.length; i++) {
		sum_Ri += R[i];
		if (G[i] > 0) {
			sum_Ri_coloc += R[i];
		}
	}
	M1 = sum_Ri_coloc / sum_Ri;
	// calculate M2
	sum_Gi_coloc = 0;
	sum_Gi = 0;
	for (i = 0; i < G.length; i++) {
		sum_Gi += G[i];
		if (R[i] > 0) {
			sum_Gi_coloc += G[i];
		}
	}
	M2 = sum_Gi_coloc / sum_Gi;
	return newArray(M1, M2);
}
