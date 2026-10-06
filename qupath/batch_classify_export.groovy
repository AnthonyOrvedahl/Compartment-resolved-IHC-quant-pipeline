/**
 * batch_classify_export.groovy
 *
 * Run: Automate > Run for project, select all images, "Save changes" on.
 * Developed with QuPath 0.7.0.
 *
 * Per image
 *   1. Print the colour-deconvolution stains currently set. They must be the
 *      cohort vectors (HEMA, DAB below). After export, confirm with
 *      pipeline/check_stains.py, which recomputes DAB OD from the pixels.
 *   2. If the image has no detections: set the image type and cohort stains,
 *      create a full-image annotation only when the image has no annotation, and
 *      run cell detection inside the annotation(s). If detections already exist,
 *      stain setting and detection are skipped.
 *   3. Classify every detection by centroid with the pixel classifier
 *      PIXEL_CLASSIFIER.
 *   4. Export every detection, all measurements, to
 *      <project>/exports/<image name>.Classifications.txt
 *
 * Notes
 *   - The DAB-positive object classifier is deliberately not run: it overwrites
 *     Classification (the tissue class) with Positive/Negative. Positivity is
 *     thresholded downstream from "Cell: DAB OD mean".
 *   - Existing annotations are never removed. When a slide carries two sections
 *     of the same block, annotate one section only; cells are detected inside the
 *     existing annotation(s).
 *   - The export reads the live measurement table, and image data is saved first.
 */

import qupath.lib.objects.classes.PathClass
import qupath.lib.gui.scripting.QPEx

def OUTDIR = buildFilePath(PROJECT_BASE_DIR, 'exports')
new File(OUTDIR).mkdirs()
def PIXEL_CLASSIFIER = 'Lung feature detection'
def HEMA = [0.71010, 0.64228, 0.28852]
def DAB  = [0.43668, 0.58566, 0.68287]

// Matches the QuPath Cell detection dialog: pixel size 0.5, expansion 7 um,
// smooth boundaries off, exclude DAB off.
def CELL_DETECTION_PARAMS = '''{
  "detectionImageBrightfield": "Hematoxylin OD",
  "requestedPixelSizeMicrons": 0.5,
  "backgroundRadiusMicrons": 8.0,
  "backgroundByReconstruction": true,
  "medianRadiusMicrons": 0.0,
  "sigmaMicrons": 1.5,
  "minAreaMicrons": 10.0,
  "maxAreaMicrons": 400.0,
  "threshold": 0.1,
  "maxBackground": 2.0,
  "watershedPostProcess": true,
  "excludeDAB": false,
  "cellExpansionMicrons": 7.0,
  "includeNuclei": true,
  "smoothBoundaries": false,
  "makeMeasurements": true
}'''

def imageData = getCurrentImageData()
def name = getProjectEntry().getImageName()
def baseName = name.replaceAll(/\.[^.]+$/, '')

// ---- 1. diagnostic: print whatever stains are currently set -----------------
def curStains = imageData.getColorDeconvolutionStains()
println name + ' | current stains: ' + (curStains ?: 'NONE SET')

// ---- 2. detect only if nothing exists yet; never remove existing annotations ----
def nCells = getDetectionObjects().size()
def existingAnns = getAnnotationObjects()
if (nCells == 0) {
    setImageType('BRIGHTFIELD_H_DAB')
    setColorDeconvolutionStains(
        '{"Name" : "cohort H-DAB", "Stain 1" : "Hematoxylin", "Values 1" : "' +
        HEMA.join(' ') + '", "Stain 2" : "DAB", "Values 2" : "' + DAB.join(' ') +
        '", "Background" : " 255 255 255"}'
    )
    if (existingAnns.isEmpty()) {
        println name + ' | no annotations found -> creating one full-image annotation (nothing to preserve here)'
        createFullImageAnnotation(true)
        existingAnns = getAnnotationObjects()
    } else {
        println name + ' | ' + existingAnns.size() + ' existing annotation(s) found - using them as-is, NOT creating a full-image annotation'
    }
    selectObjects(existingAnns)
    runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', CELL_DETECTION_PARAMS)
    getCurrentHierarchy().getSelectionModel().clearSelection()
    println name + ' | detected ' + getDetectionObjects().size() + ' cells across ' + existingAnns.size() + ' annotation(s)'
} else {
    println name + ' | ' + nCells + ' detections already present - SKIPPING stain reset and detection'
}

// ---- 3. classify every detection by centroid, cheap and idempotent ----------
def available = getProject().getPixelClassifiers().getNames()
if (!available.contains(PIXEL_CLASSIFIER)) {
    println name + ' | ERROR: pixel classifier "' + PIXEL_CLASSIFIER + '" not found. Available: ' + available
    return
}
classifyDetectionsByCentroid(PIXEL_CLASSIFIER)
fireHierarchyUpdate()

def tally = [:].withDefault { 0 }
getDetectionObjects().each { tally[it.getPathClass()?.toString() ?: 'Unclassified'] += 1 }
println name + ' | classified: ' + tally.sort { it.key }.collect { k, v -> k + ' ' + v }.join(' | ')

// ---- 4. export every detection with its full measurement set ------------------
// ObservableMeasurementTableData is the class behind Measure > Show detection
// measurements. It reads the live hierarchy; MeasurementExporter reads the saved
// .qpdata and would find nothing for unsaved detections.
def dets = getDetectionObjects()
if (dets.isEmpty()) {
    println name + ' | no detections to export'
    return
}
try {
    getProjectEntry().saveImageData(getCurrentImageData())
} catch (Exception e) {
    println name + ' | WARNING: could not save image data (' + e.getMessage() + ') - export still uses live memory'
}

def table = new qupath.lib.gui.measure.ObservableMeasurementTableData()
table.setImageData(getCurrentImageData(), dets)
def cols = table.getAllNames()
def outFile = new File(OUTDIR, baseName + '.Classifications.txt')
outFile.withWriter('UTF-8') { w ->
    w.writeLine(cols.join('\t'))
    dets.each { d ->
        w.writeLine(cols.collect { c -> table.getStringValue(d, c) ?: '' }.join('\t'))
    }
}
println name + ' | exported ' + dets.size() + ' detections, ' + cols.size() +
        ' columns -> ' + outFile.getAbsolutePath()
