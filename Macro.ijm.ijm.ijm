
inputDir = getDirectory("C:\Users\mudzimtb\Downloads\Spheroids\Day1\A");
outputCSV = inputDir + "Combined_Results.csv";

// ---------- USER PARAMETERS ----------
minSize  = 100;       // area lower bound (µm^2 if image is calibrated; else px^2)
maxSize  = 1e12;      // effectively Infinity
minCirc  = 0.00;      // 0..1
maxCirc  = 1.00;

EXCLUDE_EDGES = true; // exclude ROIs touching the image edges?
DECIMALS      = 3;    // results precision

// If your stacks are time-lapse, set frame interval (minutes). If unknown, set to 0.
FRAME_INTERVAL_MIN = 5.0;  // e.g., 5.0 minutes per frame; set 0 to skip Time column
// ------------------------------------

setBatchMode(true);
run("Clear Results");
setOption("ExpandableArrays", true);

list = getFileList(inputDir);

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".tif") && !endsWith(list[i], "_masks.tif")) {
        rawName = list[i];
        rawPath = inputDir + rawName;
        maskName = replace(rawName, ".tif", "_masks.tif");
        maskPath = inputDir + maskName;

        if (!File.exists(maskPath)) {
            print("Skipping (no mask): " + rawName);
            continue;
        }

        // --- Open raw and mask ---
        open(rawPath);
        rawTitle = getTitle();

        open(maskPath);
        maskTitle = getTitle();

        // Work on a duplicate of mask to keep original clean
        selectWindow(maskTitle);
        run("Duplicate...", "title=mask_tmp");
        close(maskTitle);
        selectWindow("mask_tmp");

        // --- Convert label mask to binary: objects (>0) = white ---
        run("32-bit");
        setThreshold(1, 1e9);           // any label >0 considered foreground
        setOption("BlackBackground", true);
        run("Convert to Mask");         // -> 8-bit binary (0/255)

        // --- Set measurements, redirect intensity to raw image ---
        selectWindow(rawTitle);
        run("Set Measurements...",
            "area mean min max centroid center perimeter shape feret redirect=["+rawTitle+"] decimal="+DECIMALS);

        // Determine if we’ve got stacks
        selectWindow("mask_tmp");
        nZ = nSlices;

        // Analyze either single slice or each slice in a stack
        if (nZ <= 1) {
            // Single image
            startRow = nResults;
            opts = "size="+minSize+"-"+maxSize+" circularity="+minCirc+"-"+maxCirc+" display add";
            if (EXCLUDE_EDGES) opts = opts + " exclude";
            run("Analyze Particles...", opts);
            endRow = nResults;

            for (r = startRow; r < endRow; r++) {
                setResult("File", r, rawName);
                setResult("Frame", r, 1);
                if (FRAME_INTERVAL_MIN > 0) setResult("Time_min", r, (1-1) * FRAME_INTERVAL_MIN);
            }
            updateResults();

        } else {
            // Stack: loop over slices and keep raw/mask in sync
            for (s = 1; s <= nZ; s++) {
                // Set slice in both images
                selectWindow("mask_tmp");
                Stack.setSlice(s);
                selectWindow(rawTitle);
                Stack.setSlice(s);

                // Analyze this slice
                startRow = nResults;
                opts = "size="+minSize+"-"+maxSize+" circularity="+minCirc+"-"+maxCirc+" display add";
                if (EXCLUDE_EDGES) opts = opts + " exclude";
                selectWindow("mask_tmp");
                run("Analyze Particles...", opts);
                endRow = nResults;

                // Annotate rows with filename + frame + time
                for (r = startRow; r < endRow; r++) {
                    setResult("File", r, rawName);
                    setResult("Frame", r, s);
                    if (FRAME_INTERVAL_MIN > 0) setResult("Time_min", r, (s-1) * FRAME_INTERVAL_MIN);
                }
                updateResults();
            }
        }

        // Cleanup windows
        close("mask_tmp");
        close(rawTitle);
    }
}

// Save one combined table
saveAs("Results", outputCSV);
setBatchMode(false);
print("Done. Saved: " + outputCSV);
