import net.imagej.Dataset;
import net.imglib2.*;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import org.scijava.command.Command;
import org.scijava.plugin.Plugin;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation>Just Cross Correlation")
public class Just_Cross_Correlation <R extends RealType<R>> extends Abstract_CCC_base {

    public Just_Cross_Correlation() {
    }

    @Override
    public void run(){
        initializePlugin(new String[]{"Original CC result"});

        //region Error checking
        dimensionChecking(dataset1, dataset2);
        multiChannelErrorCheck(dataset1,dataset2);
        //endregion

        //region Single frame analysis
        if(dataset1.getFrames() == 1) {
            try {
                radialProfiler = new RadialProfiler(dataset1, scale);
            } catch (Exception e) {
                e.printStackTrace();
                return;
            }
            try {
                colocalizationAnalysis(dataset1, dataset2, radialProfiler, intermediates, dataset1.factory());
            } catch (Exception e) {
                e.printStackTrace();
            }

            generatePlots();
            generateFullCorrelationTable();
        }
        //endregion

        //region Multi-frame analysis
        else {
            Dataset dataset1copy = dataset1.duplicate();
            Dataset dataset2copy = dataset2.duplicate();

            for (long i = 0; i < dataset1.getFrames(); i++) {
                setActiveFrame(i);
                RandomAccessibleInterval<RealType<?>> temp1 = getActiveFrame(dataset1copy);
                RandomAccessibleInterval<RealType<?>> temp2 = getActiveFrame(dataset2copy);
                if(showIntermediates){
                    for (int m = 0; m < intermediates.length; m++) {
                        intermediatesViewsPasser[m] = getActiveFrame(intermediates[m]);
                    }
                }
                try {
                    radialProfiler = new RadialProfiler(temp1, scale);
                } catch (Exception e) {
                    e.printStackTrace();
                    return;
                }

                try {
                    colocalizationAnalysis(temp1, temp2, radialProfiler, intermediatesViewsPasser, dataset1.factory());
                } catch (Exception e) {
                    e.printStackTrace();
                    throw e;
                }
                addDataToHeatmaps(i);
            }
        }

        if(showIntermediates){
            displayIntermediates();
        }
        if(saveFolder != null && saveFolder.getPath() != ""){
            saveResultsToFolder();
        }
    }

    //made this to quickly and easily test different extension methods for correlation

    private <R extends RealType<?>> void colocalizationAnalysis(RandomAccessibleInterval <R> img1, RandomAccessibleInterval <R> img2, RadialProfiler radialProfiler, RandomAccessibleInterval <R> [] localIntermediates, ImgFactory imgFactory){
        statusService.showStatus(statusBase + "Calculating original correlation");

        Img<FloatType> crossCorrelation = ops.create().img(img1, new FloatType());

        statusService.showStatus(statusBase + "Initializing data");

        CCfunctions ccFunctions = new CCfunctions(img1, img2, scale, imgFactory);

        statusService.showStatus(statusBase + "Calculating cross-correlation");

        ccFunctions.calculateCC(crossCorrelation);

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[0], crossCorrelation).multiThreaded().forEachPixel((a, b) -> a.setReal(b.get()));
        }

        statusService.showStatus(statusBase + "Calculating radial profile");
        radialProfiler.calculateOCorrProfile(crossCorrelation);
    }
}

