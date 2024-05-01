import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import org.scijava.command.Command;
import org.scijava.plugin.Plugin;
import utils.CCfunctions;
import utils.RadialProfiler;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, headless = true, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation>CCC - No confidence")
public class CCC_No_Confidence extends Abstract_CCC_gaussian {

    public CCC_No_Confidence() {
    }

    @Override
    public void run(){

        initializePlugin(new String[]{"Subtracted CC result", "Gaussian-modified CC result"});

        //region Error checking
        dimensionChecking(dataset1, dataset2);
        dimensionChecking(dataset1, maskDataset);
        multiChannelErrorCheck(dataset1,dataset2, maskDataset);
        //endregion

        //region Single frame analysis
        if(dataset1.getFrames() == 1) {
            try {
                radialProfiler = new RadialProfiler(convertedImg1, scale);
            } catch (Exception e) {
                e.printStackTrace();
                return;
            }

            try {
                colocalizationAnalysis(convertedImg1, convertedImg2, maskDataset, radialProfiler, ContributionOf1, ContributionOf2, intermediates, convertedImg1.factory());
            } catch (Exception e) {
                e.printStackTrace();
            }
            generatePlots();
            generateFullCorrelationTable();
        }
        //endregion

        //region Multi-frame analysis
        else {
            for (long i = 0; i < dataset1.getFrames(); i++) {
                setActiveFrame(i);
                RandomAccessibleInterval<FloatType> temp1 = getActiveFrame(convertedImg1);
                RandomAccessibleInterval<FloatType> temp2 = getActiveFrame(convertedImg2);
                RandomAccessibleInterval<RealType<?>> masktemp = getActiveFrame(maskDataset);
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
                    colocalizationAnalysis(temp1, temp2, masktemp, radialProfiler, ContributionOf1 == null ? null : getActiveFrame(ContributionOf1), ContributionOf2 == null ? null : getActiveFrame(ContributionOf2), intermediatesViewsPasser, convertedImg1.factory());
                } catch (Exception e) {
                    e.printStackTrace();
                    throw e;
                }
                addDataToHeatmaps(i);
                addToTimeResultsTable(i);
            }
        }
        //endregion

        generateResults();

        if(showIntermediates){
            displayIntermediates();
        }
        if(saveFolder != null && saveFolder.getPath() != ""){
            saveResultsToFolder();
        }
    }

    private <R extends RealType<?>> void colocalizationAnalysis(RandomAccessibleInterval <FloatType> img1, RandomAccessibleInterval<FloatType> img2, RandomAccessibleInterval<R> imgMask, RadialProfiler radialProfiler, final RandomAccessibleInterval <R> contribution1, final RandomAccessibleInterval <R> contribution2, RandomAccessibleInterval <R> [] localIntermediates, ImgFactory imgFactory){

        Img<FloatType> subtracted = ops.create().img(img1, new FloatType());
        Img<FloatType> gaussModifiedCorr;

        statusService.showStatus(statusBase + "Generating averaged mask");

        CCfunctions ccFunctions = new CCfunctions(img1, img2, imgMask, scale, imgFactory);

        statusService.showStatus(statusBase + "Generating subtracted correlation");

        ccFunctions.generateSubtractedCCImage(img1, img2, imgMask, subtracted, imgFactory);

        statusService.showStatus(statusBase + "Calculating radial profile");

        radialProfiler.calculateSCorrProfile(subtracted);

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[0], subtracted).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }

        if(!generateContributionImages){
            subtracted = null;
        }

        statusService.showStatus(statusBase + "Fitting gaussian to data");
        try{radialProfiler.fitGaussianCurve();}
        catch (NullPointerException e){
            generateContributionImages = false;
            logService.warn("Failed to fit gaussian curve to cross correlation of " + dataset1.getName() + " and " + dataset2.getName() + ", suggesting no correlation between the images.\nAcquired data and intermediate correlation images (if the option was selected) will still be shown. Statistical measures will be set to error values (-1).");
        }

        if(generateContributionImages) {
            statusService.showStatus(statusBase + "Determining channel contributions");
            //gaussModifiedCorr = imgFactory.create(img1);
            gaussModifiedCorr = ops.create().img(img1, new FloatType());

            ccFunctions.generateGaussianModifiedCCImage(subtracted, gaussModifiedCorr, radialProfiler);

            if(showIntermediates){
                LoopBuilder.setImages(localIntermediates[1], gaussModifiedCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
            }

            ccFunctions.calculateContributionImages(img1, img2, gaussModifiedCorr, contribution1, contribution2);
        }
    }
}

