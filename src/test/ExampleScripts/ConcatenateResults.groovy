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
import io.scif.services.DatasetIOService;
import net.imagej.axis.Axes;
import net.imagej.Dataset;
import net.imglib2.FinalInterval;
import net.imglib2.type.numeric.real.FloatType;
import org.scijava.module.ModuleInfo;
import org.scijava.table.Table;
import org.scijava.table.Tables;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;

/**
 * The purpose of this script is to concatenate CCC results from multiple images into one Table.
 * The script is NOT setup to perform most of the pre-processing steps for CCC (decon, subtract background, etc),
 * so you will likely need to perform these in advance or modify the script to work for your data.
 * The script also assumes multi-channel files with Channel 1 and 2 as the channels to evaluate
 */



#@ File[] (label="Select images to process", style="files") fileList
#@ UIService uiService
#@ OpService ops
#@ ModuleService moduleService
#@ DatasetIOService datasetioService
#@ DatasetService datasetService


Table concatenatedTable;

ArrayList<HashMap<String, Object>> concatenatedList = new ArrayList<>();
ArrayList<String> imageNames = new ArrayList<>();

ModuleInfo ccc = moduleService.getModuleById("command:CCC.Colocalization_by_Cross_Correlation");

for (int i = 0; i < fileList.length; i++) {
    Dataset originalImage = datasetioService.open(fileList[i].getPath());
    Dataset floatImage = datasetService.create(
            ops.create().imgPlus(
                    ops.convert().float32(originalImage.getImgPlus().getImg())
                    , originalImage
            )
    );

    //Use of crop here instead of hyperslice is deliberate, as the hyperslice op does not preserve metadata
    int channelAxis = floatImage.dimensionIndex(Axes.CHANNEL);
    long[] min = new long[floatImage.numDimensions()];
    long[] max = new long[floatImage.numDimensions()];
    for(int d = 0; d < min.length; ++d){
        min[d] = floatImage.min(d);
        max[d] = floatImage.max(d);
    }
    max[channelAxis] = 0;

    FinalInterval ch1Int = new FinalInterval(min, max);
    ch1 = ops.transform().crop(floatImage, ch1Int , true);

    min[channelAxis] = 1;
    max[channelAxis] = 1;
    FinalInterval ch2Int = new FinalInterval(min, max);
    ch2 = ops.transform().crop(floatImage, ch2Int , true);

    /*
    Note: for these next blocks, subtracting the median only works as a method of background subtraction if the majority
    of your pixels are background pixels, which, in my experience, is often the case with 3D images, but not often
    with 2D
     */
    ch1 = datasetService.create(
            ops.create().imgPlus(
                    ops.math().subtract(ch1,
                            new FloatType(ops.stats().median(ch1).getRealFloat())
                    ), ch1
            )
    );
    ch1 = datasetService.create(
            ops.create().imgPlus(
                    ops.math().subtract(ch2,
                            new FloatType(ops.stats().median(ch2).getRealFloat())
                    ), ch2
            )
    );


    //Edit this block to create an ideal mask for your data
    Dataset mask = datasetService.create(
            ops.create().imgPlus(
                    ops.threshold().huang(
                            ops.filter().gauss(ch1, 5)
                    ), ch1
            )
    );

    //Uncomment next three lines if you will be using the CCC autosave feature
    //ch1.setName(image.getName() + "-ch1");
    //ch2.setName(image.getName() + "-ch2");
    //mask.setName(image.getName() + "-Huang mask");

    Table tableOut = moduleService.run(ccc, false,
            "dataset1", ch1,
            "dataset2",ch2,
            "maskAbsent", false,
            "maskDataset",mask,
            "significantDigits", 4,
            "generateContributionImages",false,
            "showIntermediates",false ,
            "saveFolder", ""
    ).get().getOutput("resultsTable");

    imageNames.add(floatImage.getName());
    concatenatedList.add(
            Maps.newHashMap(
                    ImmutableMap.of("Mean", tableOut.get(0,0),
                            "StDev", tableOut.get(0,1),
                            "Confidence", tableOut.get(0,2),
                            "R-Squared", tableOut.get(0,3),
                            "Gaussian Height", tableOut.get(0,4)
                    )
            )
    );

}

concatenatedTable = Tables.wrap(concatenatedList, imageNames);

uiService.show(concatenatedTable);
