/**
 * Script to check how the Gaussian height changes based on the normalization scheme of input PSFs.
 * Attempted to see if any alterations to the normalization procedures anywhere (in pre-processing or during
 * cross-correlation) would result in the Gaussian height coming anywhere close to reflecting the count of
 * correlated particles in this over-idealized scenario. I was hoping that the normalized to sum pre-processing
 * without normalizing to the mask volume would get near this value, but this was not the case. None of the Gauss
 * height values came particularly close to 50, regardless of the normalization scheme used.
 */

import CCC.Colocalization_by_Cross_Correlation
import net.imagej.ImageJ
import org.scijava.command.CommandInfo
import org.scijava.module.ModuleInfo
import org.scijava.table.Table

ImageJ ij = new ImageJ();

ModuleInfo originalCCCinfo = new CommandInfo(Colocalization_by_Cross_Correlation.class);

ij.module().addModule(originalCCCinfo);

Correlated_Pairs correlated_pairs_original = new Correlated_Pairs(ij, "src/test/resources/green, pre-deconvolved, 10x RL.tif", "src/test/resources/green, pre-deconvolved, 10x RL.tif");
Correlated_Pairs correlated_pairs_max = new Correlated_Pairs(ij, "src/test/resources/green, normalizedToMax.tif", "src/test/resources/green, normalizedToMax.tif");
Correlated_Pairs correlated_pairs_sum = new Correlated_Pairs(ij, "src/test/resources/green, normalizedToSum.tif", "src/test/resources/green, normalizedToSum.tif");

correlated_pairs_original.generateCorrelatedImages(512, 128, 0, 50, 0);
correlated_pairs_max.generateCorrelatedImages(512, 128, 0, 50, 0);
correlated_pairs_sum.generateCorrelatedImages(512, 128, 0, 50, 0);

Table outputs_original = ij.module().run(originalCCCinfo, true, "dataset1", correlated_pairs_original.paired1, "dataset2", correlated_pairs_original.paired2, "maskAbsent", true, "maskDataset",correlated_pairs_max.paired1, "significantDigits",4, "generateContributionImages",false, "showIntermediates",false ,"saveFolder", "src/test/testOutputs/Normalization/Original").get().getOutput("resultsTable");
Table outputs_max = ij.module().run(originalCCCinfo, true, "dataset1", correlated_pairs_max.paired1, "dataset2", correlated_pairs_max.paired2, "maskAbsent", true, "maskDataset",correlated_pairs_max.paired1, "significantDigits",4, "generateContributionImages",false, "showIntermediates",false ,"saveFolder", "src/test/testOutputs/Normalization/Max").get().getOutput("resultsTable");
Table outputs_sum = ij.module().run(originalCCCinfo, true, "dataset1", correlated_pairs_sum.paired1, "dataset2", correlated_pairs_sum.paired2, "maskAbsent", true, "maskDataset",correlated_pairs_max.paired1, "significantDigits",4, "generateContributionImages",false, "showIntermediates",false ,"saveFolder", "src/test/testOutputs/Normalization/Sum").get().getOutput("resultsTable");

println("Number of correlated pairs: 50");
println("Gauss height not normalized:" + outputs_original.get(0,4));
println("Gauss height for Normalized to max:" + outputs_max.get(0,4));
println("Gauss height for Normalized to sum:" + outputs_sum.get(0,4));

ij.dispose();