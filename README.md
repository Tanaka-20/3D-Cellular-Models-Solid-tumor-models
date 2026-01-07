# Training Cellpose model for spheroids detection and segmentation.  

Spheroids are spherical cellular units that are generally cultured as free-floating aggregates and are arguably of low complexity in mirroring tumor organization. In general, organoids can be referred to as cells grown in 3D to form structural units that partially resemble the organ, both in structure and function.

## Study window: 4 days (growth observation)
Imaging: 12 plates total; per drug × day there are 12 images
- Images 1–4 = High Ca²⁺
- Images 5–8 = No spheroids formed (exclude)
- Images 9–12 = Low Ca²⁺
  Preprocessing: For each drug (A–H) and day (1–4), stacked TIFF image sequences were created.

## Purpose
Is to manually train the custom cellpose model to detect and segment spheroids and be able to have a track analysis over time be it days or minutes. 
Spheroid Formation & Growth Under Drug Targets A–H and Two Ca²⁺ Conditions.

## 1) Objective 

- Quantify how spheroid size (area) changes over 4 days across drug targets A–H under two calcium concentrations (High vs Low).
- Compare growth trajectories and infer drug × calcium effects on spheroid development.
---

## 2) Study Design & Data Map

- Factors:
	- Drug: A, B, C, D, E, F, G, H
	- Calcium: High (images 1–4), Low (images 9–12)
	- Day: 1, 2, 3, 4
- Replicates: Up to 4 images per condition (per drug × day × calcium).

  # Image Batch Processing
 ## Main Aim

To automate the segmentation and quantitative analysis of spheroid images by combining Cellpose (Python) for mask generation and Fiji/ImageJ macros for batch measurement extraction, producing one consolidated results file.

## Objectives

1. Develop a Python batch script to:
	- Apply a custom Cellpose model to all .tif images in a folder.
	- Generate segmentation masks (*_masks.tif).
	- Create preview overlays (*_masks_preview.png) for visual quality control.
2. Use a Fiji/ImageJ macro to:
	- Import both raw images and their corresponding _masks.tif.
	- Convert masks into binary ROIs suitable for analysis.
	- Run Analyze ▸ Analyze Particles with defined thresholds.
	- Redirect intensity measurements to the raw image.
	- Append metadata (Filename, Frame, Time).
	- Export one combined results table (Combined_Results.csv).
3. Streamline the entire workflow for reproducibility and scalability across multiple datasets.
##Codes Used


##Python Script (Cellpose v4 batch processing + previews) Cellpose Segmentation Script
Python script used to segment spheroidal cells using a fine-tuned Cellpose model.
The pipeline saves both raw mask outputs (`.tif`) and overlay previews (`.png`).
```python
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from cellpose import models

IMG_DIR    = Path(r"C:\Users\mudzimtb\Downloads\Spheroids\Day1\A")
PATTERN    = "*.tif"
MODEL_PATH = r"C:\Users\mudzimtb\Downloads\Spheroids\Day1\New folder\models\cpsam_20250904_165506"

DIAMETER   = None      # or e.g. 75.0
USE_GPU    = False


def save_preview(img, masks, out_path):
    """
    Save a visualization of masks over the original image.
    Grey background + colored masks overlay.
    """
    plt.figure(figsize=(6,6))
    if img.ndim == 2:
        plt.imshow(img, cmap="gray")
    elif img.ndim == 3 and img.shape[-1] in (3,4):
        plt.imshow(img)
    else:
        plt.imshow(img[0], cmap="gray")

    plt.imshow(masks, alpha=0.4, cmap="nipy_spectral")
    plt.axis("off")
    plt.savefig(out_path, dpi=200, bbox_inches="tight", pad_inches=0)
    plt.close()


def run_cellpose_v4(model, img, diameter):
    if img.ndim == 2:
        masks, *_ = model.eval(img, diameter=diameter, channels=None, channel_axis=None)
        return masks.astype(np.uint16)

    if img.ndim == 3 and img.shape[-1] in (3,4):
        masks, *_ = model.eval(img, diameter=diameter, channels=None, channel_axis=-1)
        return masks.astype(np.uint16)

    if img.ndim == 3:
        planes = []
        for k in range(img.shape[0]):
            masks, *_ = model.eval(img[k], diameter=diameter, channels=None, channel_axis=None)
            planes.append(masks.astype(np.uint16))
        return np.stack(planes, axis=0)

    raise ValueError(f"Unsupported shape {img.shape}")


def main():
    model = models.CellposeModel(pretrained_model=MODEL_PATH, gpu=USE_GPU)

    for img_path in sorted(IMG_DIR.glob(PATTERN)):
        print(f"Processing {img_path.name}...")
        img = io.imread(str(img_path))
        masks = run_cellpose_v4(model, img, DIAMETER)

        mask_out = img_path.with_name(img_path.stem + "_masks.tif")
        io.imsave(str(mask_out), masks)

        preview_out = img_path.with_name(img_path.stem + "_masks_preview.png")
        save_preview(img, masks, preview_out)

        print(f"  Saved mask: {mask_out}")
        print(f"  Saved preview: {preview_out}")


if __name__ == "__main__":
    main()
```
## Fiji / ImageJ Macro: Batch ROI Quantification

This macro processes Cellpose-generated mask images to extract
morphological and intensity measurements from corresponding raw images.
Results from all images are combined into a single CSV file.

```ijm
inputDir = getDirectory("C:\\Users\\mudzimtb\\Downloads\\Spheroids\\Day1\\A");
outputCSV = inputDir + "Combined_Results.csv";

// ---------- USER PARAMETERS ----------
minSize  = 100;       // area lower bound (µm^2 if image is calibrated; else px^2)
maxSize  = 1e12;      // effectively Infinity
minCirc  = 0.00;      // 0..1
maxCirc  = 1.00;

EXCLUDE_EDGES = true; // exclude ROIs touching the image edges?
DECIMALS      = 3;    // results precision

// If your stacks are time-lapse, set frame interval (minutes). If unknown, set to 0.
FRAME_INTERVAL_MIN = 5.0;  // minutes per frame; set 0 to skip Time column
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

        open(rawPath);
        rawTitle = getTitle();

        open(maskPath);
        maskTitle = getTitle();

        selectWindow(maskTitle);
        run("Duplicate...", "title=mask_tmp");
        close(maskTitle);
        selectWindow("mask_tmp");

        run("32-bit");
        setThreshold(1, 1e9);
        setOption("BlackBackground", true);
        run("Convert to Mask");

        selectWindow(rawTitle);
        run("Set Measurements...",
            "area mean min max centroid center perimeter shape feret redirect=["+rawTitle+"] decimal="+DECIMALS);

        selectWindow("mask_tmp");
        nZ = nSlices;

        if (nZ <= 1) {
            startRow = nResults;
            opts = "size="+minSize+"-"+maxSize+" circularity="+minCirc+"-"+maxCirc+" display add";
            if (EXCLUDE_EDGES) opts += " exclude";
            run("Analyze Particles...", opts);
            endRow = nResults;

            for (r = startRow; r < endRow; r++) {
                setResult("File", r, rawName);
                setResult("Frame", r, 1);
                if (FRAME_INTERVAL_MIN > 0)
                    setResult("Time_min", r, 0);
            }
            updateResults();

        } else {
            for (s = 1; s <= nZ; s++) {
                selectWindow("mask_tmp");
                Stack.setSlice(s);
                selectWindow(rawTitle);
                Stack.setSlice(s);

                startRow = nResults;
                opts = "size="+minSize+"-"+maxSize+" circularity="+minCirc+"-"+maxCirc+" display add";
                if (EXCLUDE_EDGES) opts += " exclude";
                selectWindow("mask_tmp");
                run("Analyze Particles...", opts);
                endRow = nResults;

                for (r = startRow; r < endRow; r++) {
                    setResult("File", r, rawName);
                    setResult("Frame", r, s);
                    if (FRAME_INTERVAL_MIN > 0)
                        setResult("Time_min", r, (s-1) * FRAME_INTERVAL_MIN);
                }
                updateResults();
            }
        }

        close("mask_tmp");
        close(rawTitle);
    }
}

saveAs("Results", outputCSV);
setBatchMode(false);
print("Done. Saved: " + outputCSV);
```


   ##Results Examples

1. For each spheroid image, the pipeline produces:
- *_masks.tif → segmentation mask (integer labels).
- *_masks_preview.png → overlay preview (raw image + colored masks).
2. Fiji macro produces Combined_Results.csv with columns:
- File → source image filename
- Area, Circularity, Perimeter, Feret’s diameter, Centroid (X,Y), Intensity metrics, etc.
  ## ROI Morphological and Intensity Measurements

| File              | Area (µm²) | Circularity | Feret Diameter | Mean Intensity |
|-------------------|------------|-------------|----------------|----------------|
| A01_x0_y0_w0.tif  | 10234.5    | 0.82        | 121.0          | 145.3          |
| A02_x0_y0_w0.tif  | 9800.7     | 0.85        | 118.9          | 150.2          |
*Quantitative morphological and fluorescence intensity measurements extracted from segmented ROIs.*

Example of a review masked image file
<img width="930" height="691" alt="image" src="https://github.com/user-attachments/assets/af687a45-595f-45d6-9432-b4eedfb1d4fa" />
This workflow provides a reproducible pipeline for going from raw images → Cellpose masks → Fiji measurements → one combined dataset ready for statistical analysis.
