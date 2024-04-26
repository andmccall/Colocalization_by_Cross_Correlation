import net.imagej.ImageJ
import org.scijava.command.CommandInfo
import org.scijava.module.ModuleInfo

ImageJ ij = new ImageJ();

ModuleInfo originalCCCinfo = new CommandInfo(Colocalization_by_Cross_Correlation.class);
ModuleInfo newCCCinfo = new CommandInfo(CCC_No_Confidence.class);

ij.module().addModule(originalCCCinfo);
ij.module().addModule(newCCCinfo);

Correlated_Pairs correlated_pairs = new Correlated_Pairs(ij);

correlated_pairs.generateCorrelatedImages(512, 128, 15, 50, 0);

originalOutputs = ij.module().run(originalCCCinfo, true, "dataset1", correlated_pairs.paired1, "dataset2", correlated_pairs.paired2, "maskAbsent", true, "maskDataset",correlated_pairs.paired1, "significantDigits",4, "generateContributionImages",false, "showIntermediates",false ,"saveFolder", "src/test/testOutputs/Accuracy");
ij.module().run(newCCCinfo, true, "dataset1", correlated_pairs.paired1, "dataset2",correlated_pairs.paired2, "maskAbsent", true, "maskDataset",correlated_pairs.paired1, "significantDigits",4, "generateContributionImages",false, "showIntermediates",false ,"saveFolder", "src/test/testOutputs/Accuracy");