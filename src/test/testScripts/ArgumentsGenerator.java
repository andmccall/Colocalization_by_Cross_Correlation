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
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import net.imagej.Dataset;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ArgumentsGenerator {

    public List<HashMap<String,Object>> jccArguments;
    public List<HashMap<String,Object>> gaussArguments;

    public ArgumentsGenerator(Dataset img1, Dataset img2, Dataset mask, File saveRoot){
        newArguments(img1, img2, mask, saveRoot);
    }
    public void newArguments(Dataset img1, Dataset img2, Dataset mask, File saveRoot){
        jccArguments = new ArrayList<>();
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", true).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", false).put("saveFolder",new File("")).build()));
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", true).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", true).put("saveFolder", new File("")).build()));
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", false).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", true).put("saveFolder", new File("")).build()));
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", true).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", false).put("saveFolder", new File(saveRoot.getAbsolutePath() + File.separator +"NoMaskNoOptions" + File.separator)).build()));
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", true).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", true).put("saveFolder", new File(saveRoot.getAbsolutePath() + File.separator +"NoMaskWithOptions" + File.separator)).build()));
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", false).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", false).put("saveFolder", new File(saveRoot.getAbsolutePath() + File.separator +"MaskNoOptions" + File.separator)).build()));
        jccArguments.add(Maps.newHashMap(ImmutableMap.<String,Object>builderWithExpectedSize(8).put("dataset1", img1).put("dataset2", img2).put("maskAbsent", false).put("maskDataset", mask).put("significantDigits", 4).put("showIntermediates", true).put("saveFolder", new File(saveRoot.getAbsolutePath() + File.separator +"MaskWithOptions" + File.separator)).build()));

        gaussArguments = new ArrayList<>();
        jccArguments.forEach((argument) ->
                gaussArguments.add((HashMap<String, Object>) argument.clone()));
        gaussArguments.get(0).put("generateContributionImages", false);
        gaussArguments.get(1).put("generateContributionImages", true);
        gaussArguments.get(2).put("generateContributionImages", true);
        gaussArguments.get(3).put("generateContributionImages", false);
        gaussArguments.get(4).put("generateContributionImages", true);
        gaussArguments.get(5).put("generateContributionImages", false);
        gaussArguments.get(6).put("generateContributionImages", true);

    }
}
