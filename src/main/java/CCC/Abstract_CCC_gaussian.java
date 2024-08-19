package CCC;

import net.imagej.Dataset;
import net.imglib2.type.numeric.real.FloatType;
import org.apache.commons.io.FileUtils;
import org.scijava.ItemIO;
import org.scijava.plugin.Parameter;
import org.scijava.table.Table;
import org.scijava.table.Tables;
import utils.RadialProfiler;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.*;

/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

//@Plugin(type = Command.class, headless = true, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation")
public abstract class Abstract_CCC_gaussian extends Abstract_CCC_base {

    @Parameter(label = "Generate contribution images?", description = "Generates images that highlight the signal from Image 1 and Image 2 that contributed to the result. Uncheck to use less memory.")
    protected boolean generateContributionImages;

    @Parameter(type = ItemIO.OUTPUT)
    protected Dataset ContributionOf1, ContributionOf2;

    @Parameter(type = ItemIO.OUTPUT, label = "CC Results")
    protected Table<org.scijava.table.Column<Double>, Double> resultsTable;

    //region Time-lapse specific variables for Gaussian fit
    @Parameter (type = ItemIO.OUTPUT, label = "Correlation over time")
    protected Table timeCorrelationTableOut;

    protected ArrayList<HashMap<String,Double>> timeCorrelationTable;
    protected ArrayList<String> timeCorrelationTableRowNames;
    //endregion

    @Override
    protected void initializePlugin(String[] intermediateNames){
        super.initializePlugin(intermediateNames);
        if (generateContributionImages) {
            initializeContributionImages();
        }
    }

    protected void initializeContributionImages(){
        ContributionOf1 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset1.getName(), inputAxisTypes);
        ContributionOf1.setAxes(inputCalibratedAxes);
        ContributionOf2 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset2.getName(), inputAxisTypes);
        ContributionOf2.setAxes(inputCalibratedAxes);
    }

    protected void generateFullCorrelationTable(){
        List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
        RadialProfiler finalRadialProfile = radialProfiler;
        SortedMap<Double, Double> keyMap;

        keyMap = radialProfiler.oCorrMap != null ? radialProfiler.oCorrMap : radialProfiler.sCorrMap;

        keyMap.keySet().stream().forEachOrdered((d) -> {
            LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
            row.put("Distance (" + getUnitType() +")", (getSigDigits(d)));
            if(radialProfiler.oCorrMap != null) row.put("Original CC", getSigDigits(finalRadialProfile.oCorrMap.get(d)));
            if(radialProfiler.sCorrMap != null) row.put("Subtracted CC", getSigDigits(finalRadialProfile.sCorrMap.get(d)));
            if(radialProfiler.gaussian != null) row.put("Gaussian fit", getSigDigits(finalRadialProfile.gaussian.value(d)));
            correlationTableList.add(row);
        });
        correlationTable = Tables.wrap(correlationTableList, null);
    }

    protected void generateResults(){
        if(dataset1.getFrames() == 1){
            generateResultsTable(radialProfiler.gaussFitParameters[1], radialProfiler.gaussFitParameters[2], radialProfiler.confidence, radialProfiler.rSquared, radialProfiler.gaussian.value(radialProfiler.gaussFitParameters[1]), null);
            addGaussianToSummaryFile(radialProfiler.gaussFitParameters[1], radialProfiler.gaussFitParameters[2], radialProfiler.confidence, radialProfiler.rSquared, radialProfiler.gaussian.value(radialProfiler.gaussFitParameters[1]), null);
        }
        else{
            HashMap<String, Double> bestFrame = getBestFrameResults();
            generateResultsTable(bestFrame.get("Mean"), bestFrame.get("SD"), bestFrame.get("Confidence"), bestFrame.get("R-squared"), bestFrame.get("Gaussian height"), timeCorrelationTable.indexOf(bestFrame));
            addGaussianToSummaryFile(bestFrame.get("Mean"), bestFrame.get("SD"), bestFrame.get("Confidence"), bestFrame.get("R-squared"), bestFrame.get("Gaussian height"), timeCorrelationTable.indexOf(bestFrame));
            timeCorrelationTableOut = Tables.wrap(timeCorrelationTable, timeCorrelationTableRowNames);
        }
    }

    protected void generateResultsTable(Double mean, Double SD, Double confidence, Double rSquared, Double gaussHeight, Integer frame){
        List<String> rowHeaders = new ArrayList<>();
        if(frame != null) rowHeaders.add("Time of best CC (" + timeCorrelationHeatMap.axis(0).unit() + ")");
        rowHeaders.add("Mean (" + getUnitType() + ")");
        rowHeaders.add("StDev (" + getUnitType() + ")");
        if(confidence != null) rowHeaders.add("Confidence");
        rowHeaders.add("R-squared");
        rowHeaders.add("Gaussian Height");

        List<Double> resultsList = new ArrayList<>();
        if(frame != null) resultsList.add(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? getSigDigits(calibratedTime.get().calibratedValue(frame)) : frame);
        resultsList.add(getSigDigits(mean));
        resultsList.add(getSigDigits(SD));
        if(confidence != null) resultsList.add(getSigDigits(confidence));
        resultsList.add(getSigDigits(rSquared));
        resultsList.add(getSigDigits(gaussHeight));

        resultsTable = Tables.wrap(resultsList, "", rowHeaders);
    }

    protected void addToTimeResultsTable(long frame) {
        if (timeCorrelationTable == null) timeCorrelationTable = new ArrayList<>();
        if (timeCorrelationTableRowNames == null) timeCorrelationTableRowNames = new ArrayList<>();

        LinkedHashMap<String, Double> gaussianMap = new LinkedHashMap<>();
        gaussianMap.put("Mean", getSigDigits(radialProfiler.gaussFitParameters[1]));
        gaussianMap.put("SD", getSigDigits(radialProfiler.gaussFitParameters[2]));
        if (radialProfiler.confidence != null) gaussianMap.put("Confidence", getSigDigits(radialProfiler.confidence));
        gaussianMap.put("R-squared", getSigDigits(radialProfiler.rSquared));
        gaussianMap.put("Gaussian height", getSigDigits(radialProfiler.gaussian.value(radialProfiler.gaussFitParameters[1])));

        timeCorrelationTable.add(gaussianMap);

        timeCorrelationTableRowNames.add((calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "" + getSigDigits(calibratedTime.get().calibratedValue(frame)) + " " + calibratedTime.get().unit() : "Frame " + frame));
    }

    protected HashMap<String, Double> getFrameResults(long frame){
        return timeCorrelationTable.get((int)frame);
    }

    protected HashMap<String, Double> getBestFrameResults(){
        Comparator<HashMap<String,Double>> tableCompare = (o1, o2) -> {
            if(o1.get("Confidence") != null){
                return o1.get("Confidence").compareTo(o2.get("Confidence"));
            }
            else{
                return o2.get("SD").compareTo(o1.get("SD"));
            }
        };
        return timeCorrelationTable.stream().max(tableCompare).isPresent() ? timeCorrelationTable.stream().max(tableCompare).get() : null;
    }



    protected void addGaussianToSummaryFile(Double mean, Double SD, Double confidence, Double rSquared, Double gaussHeight, Integer frame){
        summary = summary + "Fit a gaussian curve to the cross-correlation of: \n\""
                + dataset1.getName() +
                "\"\n with \n\"" +
                dataset2.getName() +
                "\"\n using the mask \n\"" +
                (maskAbsent? "No mask selected" : maskDataset.getName()) + "\":" +
                (frame != null ? "\n\nBest correlation found at " + (calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "time " + calibratedTime.get().calibratedValue(frame) + ", ": "") + "frame " + frame : "") +
                "\n\nMean (" + getUnitType() +"): " + getSigDigits(mean) +
                "\nStandard deviation: " + getSigDigits(SD) +
                (confidence != null ? "\n\nConfidence: " + getSigDigits(confidence) : "\n") +
                "\nR-squared: " + getSigDigits(rSquared) +
                "\n\nGaussian height (generally unused):" + getSigDigits(gaussHeight);

        if(confidence == -1){
            summary = summary + "\n\nThe negative confidence and R-squared values are an error result that indicate that a Gaussian curve could not be fit to the data.\nThis could be because no spatial correlation exists, or because the mask was too narrowly defined. The mask should NOT be a simple segmentation of the objects you want to measure.\nSee the website for details on how to generate a proper mask for your data.";
            return;
        }

        if(confidence != null && confidence < 0.1){
            summary = summary + "\n\nThe confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nSee the website for more details.";
        }

        if(rSquared != null && rSquared < 0.05){
            summary = summary + "\n\nThe R-squared value for the gaussian regression is very low.\nThis can indicate a low signal to noise ratio, or that no spatial correlation exists and the curve was fit to image noise.\nSee website for more details.";
        }
    }

    protected void saveResultsToFolder(){
        super.saveResultsToFolder();
        try {
            if(resultsTable != null) ioService.save(resultsTable,saveFolder.getAbsolutePath() + File.separator + "CC Results.csv" );
            FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + File.separator + "Summary.txt"), summary, (Charset) null);
            if(timeCorrelationTable != null){
                ioService.save(timeCorrelationTableOut, saveFolder.getAbsolutePath() + File.separator + "Gaussian fits over time.csv");
            }
            if(generateContributionImages){
                saveDatasetsToFolder(ContributionOf1, ContributionOf2);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

