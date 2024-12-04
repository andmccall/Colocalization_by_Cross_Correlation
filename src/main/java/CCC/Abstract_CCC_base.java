package CCC;

import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.LinearAxis;
import net.imagej.display.ColorTables;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;
import net.imglib2.view.Views;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.io.IOService;
import org.scijava.log.LogService;
import org.scijava.plot.*;
import org.scijava.plugin.Parameter;
import org.scijava.table.Table;
import org.scijava.table.Tables;
import org.scijava.ui.UIService;
import org.scijava.ui.swing.viewer.plot.jfreechart.XYPlotConverter;
import org.scijava.util.ColorRGB;
import utils.RadialProfiler;
import utils.ErrorChecking;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;

import static java.util.stream.Collectors.toList;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

//@Plugin(type = Command.class, headless = true, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation")
public abstract class Abstract_CCC_base implements Command{

    @Parameter
    protected LogService logService;

    @Parameter
    protected StatusService statusService;

    @Parameter
    protected UIService uiService;

    @Parameter
    protected DatasetService datasetService;

    @Parameter
    protected PlotService plotService;

    @Parameter
    protected IOService ioService;

    @Parameter
    protected DatasetIOService datasetIOService;

    @Parameter
    protected OpService ops;

    @Parameter(label = "Image 1: ", persist = false)
    protected Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    protected Dataset dataset2;

    @Parameter(label = "No mask (not recommended)?", description = "Only check this if your image has no region of interest for analysis purposes.", callback = "maskCallback")
    protected boolean maskAbsent;

    @Parameter(label = "Mask: ", description = "The mask used to define the low-spatial frequency component. This is important, more details at: imagej.net/plugins/colocalization-by-cross-correlation", required = false, persist = false)
    protected Dataset maskDataset;

    @Parameter(label = "Significant digits: ")
    protected int significantDigits;

    @Parameter(label = "Show intermediate images? ", description = "Shows images of numerous steps throughout the algorithm. Uncheck to use less memory. More details at: imagej.net/plugins/colocalization-by-cross-correlation")
    protected boolean showIntermediates;

    @Parameter(label = "Output directory (leave blank for none):", description = "The directory to automatically save all generated output, including the intermediate images if the \"Show Intermediates\" box is checked", required = false, style="directory")
    protected File saveFolder;

    @Parameter(type = ItemIO.OUTPUT)
    protected Dataset timeCorrelationHeatMap;

    @Parameter(type = ItemIO.OUTPUT, label = "Table of correlations")
    protected Table correlationTable;

    @Parameter(type = ItemIO.OUTPUT)
    protected XYPlot plot;

    protected Dataset [] intermediates;

    protected String statusBase = "";
    protected int currentStatus = 0;
    protected int maxStatus = 0;

    protected double [] scale;
    protected CalibratedAxis[] inputCalibratedAxes;
    protected AxisType[] inputAxisTypes;
    protected Img <FloatType> convertedImg1, convertedImg2;
    protected RadialProfiler radialProfiler;
    protected SCIFIOConfig config;
    protected String summary;

    //region Time-lapse variables
    protected Optional<CalibratedAxis> calibratedTime;
    protected RandomAccess<RealType<?>> correlationAccessor;
    protected long[] minDims;
    protected long[] maxDims;
    protected int timeAxis;
    protected RandomAccessibleInterval [] intermediatesViewsPasser;
    //endregion

    protected Abstract_CCC_base() {
    }

    @Override
    protected void finalize() throws Throwable {
        statusService.showStatus(maxStatus, maxStatus,statusBase + "Finished!");
        super.finalize();
    }

    protected void initializePlugin(String[] intermediateNames){
        statusService.showStatus("Initializing plugin data");

        //region Plugin initialization (mostly creating datasets)
        String version = "2.2.3";
        summary = "Results generated using CCC version " + version + "\n"
                + "Plugin website: https://imagej.net/plugins/colocalization-by-cross-correlation\n\n";

        if(saveFolder != null) saveFolder.mkdirs();

        config = new SCIFIOConfig();
        config.writerSetFailIfOverwriting(false);

        inputCalibratedAxes = new CalibratedAxis[dataset1.numDimensions()];
        inputAxisTypes = new AxisType[dataset1.numDimensions()];
        for (int i = 0; i < dataset1.numDimensions(); ++i) {
            inputCalibratedAxes[i] = dataset1.axis(i);
            inputAxisTypes[i] = dataset1.axis(i).type();
        }
        if (showIntermediates) {
            intermediates = new Dataset[intermediateNames.length];
            for (int i = 0; i < intermediateNames.length; i++) {
                intermediates[i] = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), intermediateNames[i], inputAxisTypes);
                intermediates[i].setAxes(inputCalibratedAxes);
            }
        }

        convertedImg1 = ops.convert().float32((Img) dataset1.getImgPlus());
        convertedImg2 = ops.convert().float32((Img) dataset2.getImgPlus());

        //Zero all the data outside the image mask, to prevent it from contributing to the cross-correlation result.
        if(!maskAbsent) {
            statusService.showStatus(statusBase + "Applying masks");
            LoopBuilder.setImages(convertedImg1, maskDataset).multiThreaded().forEachPixel((a, b) -> {
                if ((b.getRealDouble() == 0.0)) {
                    a.setReal(b.getRealDouble());
                }
            });
            LoopBuilder.setImages(convertedImg2, maskDataset).multiThreaded().forEachPixel((a, b) -> {
                if ((b.getRealDouble() == 0.0)) {
                    a.setReal(b.getRealDouble());
                }
            });
        }
        else{
            maskDataset = dataset1.duplicateBlank();
            LoopBuilder.setImages(maskDataset).multiThreaded().forEachPixel(SetOne::setOne);
        }

        //region Error checking
        ErrorChecking.dimensionChecking(dataset1, dataset2);
        ErrorChecking.dimensionChecking(dataset1, maskDataset);
        ErrorChecking.multiChannelErrorCheck(dataset1,dataset2, maskDataset);
        //endregion

        if(dataset1.getFrames() == 1){
            scale = new double[dataset1.numDimensions()];
            for (int i = 0; i < scale.length; i++) {
                scale[i] = dataset1.averageScale(i);
            }
        }
        else {
            maxStatus = maxStatus * (int)dataset1.getFrames();
            timeAxis = dataset1.dimensionIndex(Axes.TIME);
            calibratedTime = dataset1.axis(Axes.TIME);
            minDims = new long[dataset1.numDimensions()];
            maxDims = dataset1.dimensionsAsLongArray();

            scale = new double[dataset1.numDimensions()-1];
            int j = 0;
            for (int i = 0; i < dataset1.numDimensions(); i++) {
                if(i != timeAxis) {
                    maxDims[i] = maxDims[i]-1;
                    scale[j++] = dataset1.averageScale(i);
                }
            }
            if(showIntermediates){
                intermediatesViewsPasser = new RandomAccessibleInterval[intermediates.length];
            }
        }
    }

    protected String getUnitType(){ return dataset1.axis(Axes.X).isPresent() ? dataset1.axis(Axes.X).get().unit(): "Unlabeled distance unit"; }

    protected double getSigDigits(double input){
        BigDecimal bd = new BigDecimal(input);
        bd = bd.round(new MathContext(significantDigits));
        return bd.doubleValue();
    }

    protected void setActiveFrame(long frame){
        statusBase = "Frame " + (frame+1) + " - ";

        minDims[timeAxis] = frame;
        maxDims[timeAxis] = frame;
    }

    protected RandomAccessibleInterval<FloatType> getActiveFrame(Img<FloatType> in){
        return Views.dropSingletonDimensions(Views.interval(in, minDims, maxDims));
    }

    protected RandomAccessibleInterval<RealType<?>> getActiveFrame(Dataset in){
        return Views.dropSingletonDimensions(Views.interval(in, minDims, maxDims));
    }

    protected void generatePlots(){

        SeriesStyle oCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("blue"), LineStyle.NONE, MarkerStyle.CIRCLE);
        SeriesStyle sCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("green"), LineStyle.NONE, MarkerStyle.CIRCLE);
        SeriesStyle gaussStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("#ff00ff"), LineStyle.NONE, MarkerStyle.FILLEDCIRCLE);

        plot = plotService.newXYPlot();

        plot.xAxis().setLabel("Distance (" + getUnitType() + ")");
        plot.yAxis().setLabel("Cross correlation");
        plot.setTitle("Correlation of images");

        if(radialProfiler.gaussian != null) {
            XYSeries gaussData = plot.addXYSeries();
            gaussData.setValues(new ArrayList<>(radialProfiler.sCorrMap.keySet()), radialProfiler.sCorrMap.keySet().stream().map(radialProfiler.gaussian::value).collect(toList()));
            gaussData.setStyle(gaussStyle);
            gaussData.setLabel("Gaussian Fit");

            plot.xAxis().setManualRange(0, Math.max(radialProfiler.gaussFitParameters[1] + (5 * radialProfiler.gaussFitParameters[2]), 0.01));
        }

        if(radialProfiler.sCorrMap != null) {
            XYSeries sCorrPlotData = plot.addXYSeries();
            sCorrPlotData.setValues(new ArrayList<>(radialProfiler.sCorrMap.keySet()), new ArrayList<>(radialProfiler.sCorrMap.values()));
            sCorrPlotData.setStyle(sCorrStyle);
            sCorrPlotData.setLabel("Subtracted CC");
        }

        if(radialProfiler.oCorrMap != null) {
            XYSeries oCorrPlotData = plot.addXYSeries();
            oCorrPlotData.setValues(new ArrayList<>(radialProfiler.oCorrMap.keySet()), new ArrayList<>(radialProfiler.oCorrMap.values()));
            oCorrPlotData.setStyle(oCorrStyle);
            oCorrPlotData.setLabel("Original CC");
        }
    }

    protected void generateFullCorrelationTable(){
        List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
        RadialProfiler finalRadialProfile = radialProfiler;
        SortedMap<Double, Double> keyMap = radialProfiler.oCorrMap != null ? radialProfiler.oCorrMap : radialProfiler.sCorrMap;

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

    protected void addDataToHeatmaps(long frame){
        SortedMap<Double, Double> keyMap = radialProfiler.oCorrMap != null ? radialProfiler.oCorrMap : radialProfiler.sCorrMap;
        if(timeCorrelationHeatMap == null) {
            int channelCount = 0;
            if(radialProfiler.oCorrMap != null) ++channelCount;
            if(radialProfiler.sCorrMap != null) ++channelCount;
            if(radialProfiler.gaussian != null) ++channelCount;

            timeCorrelationHeatMap = datasetService.create(new FloatType(), new long[]{dataset1.dimension(Axes.TIME), keyMap.keySet().size(), channelCount}, "Correlation over time of " + dataset1.getName() + " and " + dataset2.getName(), new AxisType[]{Axes.X, Axes.Y, Axes.CHANNEL});

            ((LinearAxis) timeCorrelationHeatMap.axis(0)).setScale(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? calibratedTime.get().calibratedValue(1) : 1);
            timeCorrelationHeatMap.axis(0).setUnit((dataset1.axis(Axes.TIME).isPresent() ? dataset1.axis(Axes.TIME).get().unit() : "frame"));

            ((LinearAxis) timeCorrelationHeatMap.axis(1)).setScale((keyMap.lastKey()- keyMap.firstKey())/ keyMap.keySet().size());
            ((LinearAxis) timeCorrelationHeatMap.axis(1)).setOrigin(keyMap.firstKey());
            timeCorrelationHeatMap.axis(1).setUnit(dataset1.axis(Axes.X).get().unit());

            timeCorrelationHeatMap.axis(2).setType(Axes.CHANNEL);
            timeCorrelationHeatMap.initializeColorTables(channelCount);
            timeCorrelationHeatMap.setColorTable(ColorTables.MAGENTA,0);
            //EnumeratedAxis seems to be broken and doesn't show the proper values in ImageJ
            //timeCorrelationHeatMap.setAxis(new EnumeratedAxis(Axes.Y, getUnitType(), radialProfile.oCorrMap.keySet().stream().mapToDouble(Double::doubleValue).toArray()), 1);
            //uiService.showDialog("First key: " + radialProfile.oCorrMap.firstKey() + "," + timeCorrelationHeatMap.axis(1).rawValue(radialProfile.oCorrMap.firstKey()));

            correlationAccessor = timeCorrelationHeatMap.randomAccess();
        }

        double[] keySet = keyMap.keySet().stream().mapToDouble(Double::doubleValue).toArray();

        for (int k = 0; k < keySet.length; k++) {
            int currentChannel = 0;
            if(radialProfiler.gaussian != null){
                correlationAccessor.setPosition(new long[]{frame,k,currentChannel++});
                correlationAccessor.get().setReal(radialProfiler.gaussian.value(keySet[k]));
            }
            if(radialProfiler.sCorrMap != null){
                correlationAccessor.setPosition(new long[]{frame,k,currentChannel++});
                correlationAccessor.get().setReal(radialProfiler.sCorrMap.get(keySet[k]));
            }
            if(radialProfiler.oCorrMap != null){
                correlationAccessor.setPosition(new long[]{frame,k,currentChannel++});
                correlationAccessor.get().setReal(radialProfiler.oCorrMap.get(keySet[k]));
            }
        }
    }

    protected void displayIntermediates(){
        for (Dataset intermediate : intermediates) {
            uiService.show(intermediate);
        }
    }

    protected void saveDatasetsToFolder(Dataset... input){
        for (Dataset toSave: input) {
            try {
                datasetIOService.save(toSave, saveFolder.getAbsolutePath() + File.separator + toSave.getName() + ".tif", config);
            }
            catch (IOException e){
                logService.error("Exception during dataset saving.");
                logService.error(e);
            }
        }
    }

    protected void saveResultsToFolder(){
        if(!saveFolder.exists() || !saveFolder.canWrite()){
            logService.error("Output directory does not exist or does not have write permissions");
            return;
        }
        if(dataset1.getFrames() ==1) {
            try {
                File plotout = new File(saveFolder.getAbsolutePath() + File.separator + plot.getTitle() + ".png");
                XYPlotConverter converter = new XYPlotConverter();
                ChartUtils.saveChartAsPNG(plotout, converter.convert(plot, JFreeChart.class), plot.getPreferredWidth() * 2, plot.getPreferredHeight() * 2);

                ioService.save(correlationTable, saveFolder.getAbsolutePath() + File.separator + plot.getTitle() + ".csv");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        else{
            saveDatasetsToFolder(timeCorrelationHeatMap);
        }
        if (showIntermediates) {
            saveDatasetsToFolder(intermediates);
        }
    }
}

