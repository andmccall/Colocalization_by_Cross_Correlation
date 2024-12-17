/*-
 * #%L
 * Scijava plugin for spatial correlation
 * %%
 * Copyright (C) 2019 - 2024 Andrew McCall, University at Buffalo
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package CCC;

import net.imglib2.*;
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

@Plugin(type = Command.class, headless = true, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation>CCC")
public class Colocalization_by_Cross_Correlation extends Abstract_CCC_gaussian {

    public Colocalization_by_Cross_Correlation() {
    }

    @Override
    public void run(){
        maxStatus = 5;
        initializePlugin(new String[]{"Original CC result", "Subtracted CC result", "Gaussian-modified CC result"});

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
        if(saveFolder != null && !saveFolder.getPath().equals("")){
            saveResultsToFolder();
        }
    }

    private <R extends RealType<?>> void colocalizationAnalysis(RandomAccessibleInterval <FloatType> img1, RandomAccessibleInterval<FloatType> img2, RandomAccessibleInterval<R> imgMask, RadialProfiler radialProfiler, final RandomAccessibleInterval <R> contribution1, final RandomAccessibleInterval <R> contribution2, RandomAccessibleInterval <R> [] localIntermediates, ImgFactory<FloatType> floatTypeImgFactory){
        Img<FloatType> oCorr = ops.create().img(img1, new FloatType());
        Img<FloatType> subtracted = ops.create().img(img1, new FloatType());
        Img<FloatType> gaussModifiedCorr;

        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Generating averaged mask");

        CCfunctions ccFunctions = new CCfunctions(img1, img2, imgMask, scale, floatTypeImgFactory);

        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Calculating original correlation");

        ccFunctions.calculateCC(oCorr);

        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Generating subtracted correlation");

        ccFunctions.generateSubtractedCCImage(img1, img2, imgMask, subtracted, floatTypeImgFactory);

        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Calculating radial profile");
        radialProfiler.calculateBothProfiles(oCorr, subtracted);

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[0], oCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }
        oCorr = null;

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[1], subtracted).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }

        if(!generateContributionImages){
            subtracted = null;
        }

        statusService.showStatus(currentStatus++, maxStatus,statusBase + "Fitting gaussian to data");
        try{radialProfiler.fitGaussianCurve();}
        catch (NullPointerException e){
            logService.warn("Failed to fit gaussian curve to cross correlation of " + dataset1.getName() + " and " + dataset2.getName() + ", suggesting no correlation between the images.\nAcquired data and intermediate correlation images (if the option was selected) will still be shown. Statistical measures will be set to error values (-1).");
            if(generateContributionImages)
                ++currentStatus;
            return;
        }

        if(generateContributionImages) {
            statusService.showStatus(currentStatus++, maxStatus,statusBase + "Determining channel contributions");
            //gaussModifiedCorr = imgFactory.create(img1);
            gaussModifiedCorr = ops.create().img(img1, new FloatType());

            ccFunctions.generateGaussianModifiedCCImage(subtracted, gaussModifiedCorr, radialProfiler.gaussian, radialProfiler.gaussFitParameters);

            if(showIntermediates){
                LoopBuilder.setImages(localIntermediates[2], gaussModifiedCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
            }

            ccFunctions.calculateContributionImages(img1, img2, gaussModifiedCorr, contribution1, contribution2);
        }
    }
}

