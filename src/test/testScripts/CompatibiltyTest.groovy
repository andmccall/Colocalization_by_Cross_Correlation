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
import CCC.CCC_No_Confidence
import CCC.Colocalization_by_Cross_Correlation
import CCC.Just_Cross_Correlation
import net.imagej.Dataset
import net.imagej.ImageJ
import net.imagej.axis.Axes
import net.imglib2.view.Views
import org.scijava.command.CommandInfo
import org.scijava.module.ModuleInfo

ImageJ ij = new ImageJ();

ModuleInfo [] commands = new ModuleInfo[3];
//value of true shows Scijava outputs, value of false isn't supposed to, is for some unknown reason;
boolean scijavaProcess = false;

commands[0] = new CommandInfo(Colocalization_by_Cross_Correlation.class);
commands[1] = new CommandInfo(CCC_No_Confidence.class);
commands[2] = new CommandInfo(Just_Cross_Correlation.class);

commands.each {command ->
    ij.module().addModule(command);
}

//region: 2D 8-bit image test
Dataset flatImage = ij.scifio().datasetIO().open("http://imagej.net/images/FluorescentCells.jpg");

int channelAxis = flatImage.dimensionIndex(Axes.CHANNEL);
Dataset ch1 = ij.dataset().create(Views.hyperSlice(flatImage,channelAxis, 0));
Dataset ch2 = ij.dataset().create(Views.hyperSlice(flatImage,channelAxis, 1));
Dataset mask = ij.dataset().create(ij.op().threshold().huang(ch1));

ch1.setName("ch1");
ch2.setName("ch2");
mask.setName("mask");

ArgumentsGenerator argumentsGenerator = new ArgumentsGenerator(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/2D_8-bit/JCC"));
argumentsGenerator.jccArguments.each{ argumentArray ->
    ij.module().run(commands[2], scijavaProcess, argumentArray).get(); //note: using .get() here forces the script to wait for the results. Running multiple commands in parallel can cause memory issues
}
argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/2D_8-bit/CCC"));
argumentsGenerator.gaussArguments.each { argumentArray ->
    ij.module().run(commands[0], scijavaProcess, argumentArray).get();
}
argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/2D_8-bit/NoCon"));
argumentsGenerator.gaussArguments.each { argumentArray ->
    ij.module().run(commands[1], scijavaProcess, argumentArray).get();
}
//endregion

//region: 3D 16-bit image test
Dataset corti = ij.scifio().datasetIO().open("src/test/resources/organ-of-corti.tif");
channelAxis = corti.dimensionIndex(Axes.CHANNEL);
ch1 = ij.dataset().create(Views.hyperSlice(corti,channelAxis, 0));
ch2 = ij.dataset().create(Views.hyperSlice(corti,channelAxis, 1));
mask = ij.dataset().create(ij.op().threshold().huang(ij.op().filter().gauss(ch1, 10)));

argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/3D_16-bit/JCC"));
argumentsGenerator.jccArguments.each { argumentArray ->
    ij.module().run(commands[2], scijavaProcess, argumentArray).get();
}
argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/3D_16-bit/CCC"));
argumentsGenerator.gaussArguments.each { argumentArray ->
    ij.module().run(commands[0], scijavaProcess, argumentArray).get();
}
argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/3D_16-bit/NoCon"));
argumentsGenerator.gaussArguments.each { argumentArray ->
    ij.module().run(commands[1], scijavaProcess, argumentArray).get();
}

//endregion

//region: 4D 16-bit test
Dataset mitosis = ij.scifio().datasetIO().open("http://imagej.net/images/mitosis.tif");
channelAxis = mitosis.dimensionIndex(Axes.CHANNEL);
ch1 = ij.dataset().create(Views.hyperSlice(mitosis,channelAxis, 0));
ch2 = ij.dataset().create(Views.hyperSlice(mitosis,channelAxis, 1));
mask = ij.dataset().create(ij.op().threshold().huang(ij.op().filter().gauss(ch1, 5)));

    argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/4D_16-bit/JCC"));
    argumentsGenerator.jccArguments.each { argumentArray ->
        ij.module().run(commands[2], scijavaProcess, argumentArray).get();
    }
    argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/4D_16-bit/CCC"));
    argumentsGenerator.gaussArguments.each { argumentArray ->
        ij.module().run(commands[0], scijavaProcess, argumentArray).get();
    }
    argumentsGenerator.newArguments(ch1, ch2, mask, new File("src/test/testOutputs/Compatibility/4D_16-bit/NoCon"));
    argumentsGenerator.gaussArguments.each { argumentArray ->
        ij.module().run(commands[1], scijavaProcess, argumentArray).get();
    }
//endregion

//todo: need to add 32-bit image test for completeness.

ij.dispose();
return;
