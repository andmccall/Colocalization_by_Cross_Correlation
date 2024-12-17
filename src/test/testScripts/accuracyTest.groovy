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
import CCC.Colocalization_by_Cross_Correlation
import net.imagej.ImageJ
import org.scijava.command.CommandInfo
import org.scijava.module.ModuleInfo
import org.scijava.table.Table

ImageJ ij = new ImageJ();

ModuleInfo originalCCCinfo = new CommandInfo(Colocalization_by_Cross_Correlation.class);

ij.module().addModule(originalCCCinfo);

Correlated_Pairs correlated_pairs = new Correlated_Pairs(ij, "src/test/resources/green, pre-deconvolved, 10x RL.tif", "src/test/resources/red, pre-deconvolved, 10x RL.tif");

correlated_pairs.generateCorrelatedImages(512, 128, 15, 50, 0);


Table outputs = ij.module().run(originalCCCinfo, true, "dataset1", correlated_pairs.paired1, "dataset2", correlated_pairs.paired2, "maskAbsent", true, "maskDataset",correlated_pairs.paired1, "significantDigits",4, "generateContributionImages",false, "showIntermediates",false ,"saveFolder", "src/test/testOutputs/Accuracy").get().getOutput("resultsTable");
println((outputs.get(0,0)/15) *100 + "% accuracy\n");

ij.dispose();
