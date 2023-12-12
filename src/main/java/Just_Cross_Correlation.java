import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.LinearAxis;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccess;
import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.math.ImgMath;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import org.apache.commons.io.FileUtils;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.io.IOService;
import org.scijava.log.LogService;
import org.scijava.plot.*;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.table.Table;
import org.scijava.table.Tables;
import org.scijava.ui.DialogPrompt;
import org.scijava.ui.UIService;
import org.scijava.ui.swing.viewer.plot.jfreechart.XYPlotConverter;
import org.scijava.util.ColorRGB;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.charset.Charset;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, menuPath = "Analyze>Colocalization>Just Cross Correlation")
public class Just_Cross_Correlation implements Command{

    @Parameter
    private LogService logService;

    @Parameter
    private StatusService statusService;

    @Parameter
    private UIService uiService;

    @Parameter
    private DatasetService datasetService;

    @Parameter
    private PlotService plotService;

    @Parameter
    private IOService ioService;

    @Parameter
    private DatasetIOService datasetIOService;

    @Parameter
    private OpService ops;

    @Parameter(label = "Image 1: ", persist = false)
    private Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    private Dataset dataset2;

    @Parameter(label = "Significant digits: ")
    private int significantDigits;

    @Parameter(label = "Output directory (leave blank for none):", description = "The directory to automatically save all generated output, including the intermediate images if the \"Show Intermediates\" box is checked", required = false, style="directory")
    private File saveFolder;

    @Parameter(type = ItemIO.OUTPUT)
    private Dataset timeCorrelationHeatMap;

    @Parameter(type = ItemIO.OUTPUT, label = "Table of correlations")
    private Table correlationTable;

    @Parameter(type = ItemIO.OUTPUT, label = "CC Results")
    private Table results;

    @Parameter(type = ItemIO.OUTPUT)
    private XYPlot plot;

    @Parameter(type = ItemIO.OUTPUT, label = "Summary of Results")
    private String notes;

    private double sigDigits;

    public SortedMap<Double, Double> CorrMap;

    private double [] scale;

    private String statusBase = "";

    public Just_Cross_Correlation() {
    }

    @Override
    public void run(){

        //region Error checking
        if(dataset1.numDimensions() != dataset2.numDimensions() || dataset1.getFrames() != dataset2.getFrames()){
            logService.error("All image dimensions (XYZ, and time) must match");
            return;
        }

        if(dataset1.getChannels() > 1 || dataset2.getChannels() > 1){
            logService.error("Multi-channel images are not supported, requires separate channels");
            return;
        }
        //endregion

        statusService.showStatus("Initializing plugin data");

        //region Plugin initialization (mostly creating datasets)
        RadialProfiler radialProfile;
        CorrMap = new TreeMap<>();
        sigDigits = Math.pow(10.0, significantDigits);

        SCIFIOConfig config = new SCIFIOConfig();
        config.writerSetFailIfOverwriting(false);

        // Cannot use duplicateBlank() for creating the upcoming images, as they need to be 32-bit Float images
        CalibratedAxis [] calibratedAxes = new CalibratedAxis[dataset1.numDimensions()];
        AxisType [] axisTypes = new AxisType[dataset1.numDimensions()];
        for (int i = 0; i < dataset1.numDimensions(); ++i) {
            calibratedAxes[i] = dataset1.axis(i);
            axisTypes[i] = dataset1.axis(i).type();
        }

        //endregion

        //region Single frame analysis
        if(dataset1.getFrames() == 1) {
            scale = new double[dataset1.numDimensions()];
            for (int i = 0; i < scale.length; i++) {
                scale[i] = dataset1.averageScale(i);
            }

            try {
                radialProfile = new RadialProfiler(dataset1, scale);
            } catch (Exception e) {
                e.printStackTrace();
                return;
            }
            try{colocalizationAnalysis(dataset1.duplicate(), dataset2.duplicate(),radialProfile);}
            catch (Exception e){
                e.printStackTrace();
                throw e;
            }

            SeriesStyle oCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("blue"), LineStyle.SOLID, MarkerStyle.NONE);

            plot = plotService.newXYPlot();

            XYSeries oCorrPlotData = plot.addXYSeries();
            oCorrPlotData.setValues(CorrMap.keySet().stream().mapToDouble(Double::doubleValue).boxed().collect(Collectors.toList()), new ArrayList<>(CorrMap.values()));
            oCorrPlotData.setStyle(oCorrStyle);
            oCorrPlotData.setLabel("Original CC");

            plot.xAxis().setLabel("Distance (" + getUnitType() + ")");

            plot.yAxis().setLabel("Cross correlation");

            plot.setTitle("Correlation of images");


            List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
            CorrMap.keySet().stream().forEachOrdered((d) -> {
                LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
                row.put("Distance (" + getUnitType() +")", (getSigDigits(d)));
                row.put("Cross-correlation", getSigDigits(CorrMap.get(d)));
                correlationTableList.add(row);
            });
            correlationTable = Tables.wrap(correlationTableList, null);

            if(saveFolder != null){
                if(!saveFolder.exists() || !saveFolder.canWrite()){
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);

                    File plotout = new File(saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".png");
                    XYPlotConverter converter = new XYPlotConverter();
                    ChartUtils.saveChartAsPNG(plotout, converter.convert(plot, JFreeChart.class), plot.getPreferredWidth()*2, plot.getPreferredHeight()*2);

                    ioService.save(correlationTable, saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".csv");


                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        //endregion

        //region Multi-frame analysis
        else{
            long[] min = new long[dataset1.numDimensions()];
            long[] max = dataset1.dimensionsAsLongArray();
            int timeAxis = dataset1.dimensionIndex(Axes.TIME);

            Optional<CalibratedAxis> calibratedTime = dataset1.axis(Axes.TIME);

            scale = new double[dataset1.numDimensions()-1];
            int j = 0;
            for (int i = 0; i < dataset1.numDimensions(); i++) {
                if(i != timeAxis) {
                    scale[j++] = dataset1.averageScale(i);
                }
            }

            for (int i = 0; i < max.length; i++) {
                max[i] = max[i]-1;
            }

            min[timeAxis] = 0;
            max[timeAxis] = 0;

            RandomAccess<RealType<?>> correlationAccessor = null;

            Dataset tempHeatMap = null;

            for (long i = 0; i < dataset1.getFrames(); i++) {
                statusService.showProgress((int)i, (int)dataset1.getFrames());
                //StatusService message update is performed in colocalizationAnalysis method
                statusBase = "Frame " + (i+1) + " - ";
                min[timeAxis] = i;
                max[timeAxis] = i;

                RandomAccessibleInterval temp1 = Views.dropSingletonDimensions(Views.interval(dataset1, min, max));
                RandomAccessibleInterval temp2 = Views.dropSingletonDimensions(Views.interval(dataset2, min, max));

                try {
                    radialProfile = new RadialProfiler(temp1, scale);
                } catch (Exception e) {
                    e.printStackTrace();
                    return;
                }

                try{colocalizationAnalysis(datasetService.create(temp1), datasetService.create(temp2), radialProfile);}
                catch (Exception e){
                    e.printStackTrace();
                    throw e;
                }

                if(tempHeatMap == null) {
                    tempHeatMap = datasetService.create(new FloatType(), new long[]{dataset1.dimension(Axes.TIME), CorrMap.keySet().size()}, "Correlation over time of " + dataset1.getName() + " and " + dataset2.getName(), new AxisType[]{Axes.X, Axes.Y});
                    correlationAccessor = tempHeatMap.randomAccess();

                    ((LinearAxis)tempHeatMap.axis(1)).setScale(CorrMap.keySet().stream().mapToDouble(Double::doubleValue).max().getAsDouble()/CorrMap.keySet().size());
                }

                double[] keySet = CorrMap.keySet().stream().mapToDouble(Double::doubleValue).toArray();

                for (int k = 0; k < keySet.length; k++) {
                    correlationAccessor.setPosition(new long[]{i,k});
                    correlationAccessor.get().setReal(CorrMap.get(radialProfile.getBD(keySet[k]).doubleValue()));
                }

            }

            timeCorrelationHeatMap = tempHeatMap.copy();

            ((LinearAxis) timeCorrelationHeatMap.axis(0)).setScale(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? calibratedTime.get().calibratedValue(1) : 1);
            timeCorrelationHeatMap.setAxis(tempHeatMap.axis(1).copy(), 1);

            timeCorrelationHeatMap.axis(0).setUnit((dataset1.axis(Axes.TIME).isPresent() ? dataset1.axis(Axes.TIME).get().unit() : "frame"));
            timeCorrelationHeatMap.axis(0).setType(Axes.X);
            timeCorrelationHeatMap.axis(1).setUnit(getUnitType());
            timeCorrelationHeatMap.axis(1).setType(Axes.Y);

            timeCorrelationHeatMap.setName("Heat map of correlation over time between " + dataset1.getName() + " and " + dataset2.getName());

            if(saveFolder != null) {
                if (!saveFolder.exists() || !saveFolder.canWrite()) {
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);

                    datasetIOService.save(timeCorrelationHeatMap, saveFolder.getAbsolutePath() + "\\" + timeCorrelationHeatMap.getName(), config);

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        //endregion
    }

    private String getUnitType(){
        return dataset1.axis(Axes.X).isPresent() ? dataset1.axis(Axes.X).get().unit(): "Unlabeled distance unit";
    }

    //made this to quickly and easily test different extension methods for correlation

    private <T extends RealType> void colocalizationAnalysis(Img <T> img1, Img <T> img2, RadialProfiler radialProfiler){
        statusService.showStatus(statusBase + "Calculating original correlation");

        ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
        Img<FloatType> crossCorrelation = imgFactory.create(img1);

        statusService.showStatus(statusBase + "Initializing data");

        CCfunctions ccFunctions = new CCfunctions(img1, img2, scale);

        statusService.showStatus(statusBase + "Calculating cross-correlation");

        ccFunctions.calculateCC(crossCorrelation);

        statusService.showStatus(statusBase + "Calculating radial profile");
        radialProfiler.calculateSingleProfile(crossCorrelation, CorrMap);
    }

    private double getSigDigits(double input){ return ((Math.round(input* sigDigits))/ sigDigits);}


    //troubleshooting method for showing images at key points
/*    private void showScaledImg(Img input, String title){
        uiService.show(title, input);
    }*/

}

