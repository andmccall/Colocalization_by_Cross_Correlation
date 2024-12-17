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
package utils;

import net.imagej.Dataset;

import java.util.InputMismatchException;

public class ErrorChecking {
    public static void dimensionChecking(Dataset check1, Dataset check2){
        if(check1.numDimensions() != check2.numDimensions() || check1.getHeight() != check2.getHeight() || check1.getWidth() != check2.getWidth() || check1.getDepth() != check2.getDepth() || check1.getFrames() != check2.getFrames()){
            throw new InputMismatchException("Dimensions of " + check1.getName() + " and " + check2.getName() + " do not match.");
        }
    }

    public static void multiChannelErrorCheck(Dataset... setOfChecks){
        for (Dataset check: setOfChecks){
            if(check.getChannels() > 1){
                throw new IllegalArgumentException("Multi-channel images are not supported, requires separate channels");
            }
        }
    }
}
