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
