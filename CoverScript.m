%% FBA-based gene selection
% Please enter "no" if this isn't the first attempt
runFBA=input('Would you like to perform FBA-based gene selection?\nPlease input "yes" or "no".\n','s');
if strcmp(runFBA,'yes')
    disp('Please name the file containing your metabolic model "model.mat".')
    FBA;
end

%% This section may be skipped if this isn't the first attempt
disp('If you did not perform FBA-based gene selection, the genes that will subsequently be displayed are from the default model E. coli iAF1260.');
load FBAresult.mat;
disp(FBAresult)
Input;

%% Configurations
%The number of cells of each type in this configuration may be altered by
%changing the values assigned to the variables cellsoftype1 and
%cellsoftype2 in the script "Random.m".
%The default number of iterations for each script is 100, corresponding to
%1000 minutes.
runRandom=input('Would you like to run a configuration with cells distributed randomly? Input "yes" or "no".\n','s');
if strcmp(runRandom,'yes')
    Random;
    disp('The result has been saved in the video file "Random.avi".')
end

runLine=input('Would you like to run a configuration with cells in parallel lines? Input "yes" or "no".\n','s');
if strcmp(runLine,'yes')
    Lines;
    disp('The result has been saved in the video file "Lines.avi".')
end

runSingleLine=input('Would you like to run a configuration with cells of two types alternating in a single line? Input "yes" or "no".\n','s');
if strcmp(runSingleLine,'yes')
    SingleLine;
    disp('The result has been saved in the video file "SingleLine.avi".')
end

runConcentric=input('Would you like to run a configuration with cells of two types in concentric circles? Input "yes" or "no".\n','s');
if strcmp(runConcentric,'yes')
    Concentric;
    disp('The result has been saved in the video file "Concentric.avi".')
end

runSquare=input('Would you like to run a configuration with cells of the same type at opposite corners of a square? Input "yes" or "no".\n','s');
if strcmp(runSquare,'yes')
    Square;
    disp('The result has been saved in the video file "Square.avi".')
end

%A script without any cell locations having been assigned has been provided
%titled "Blank.m". You may place cells in any configuration of interest in
%this script and run it.