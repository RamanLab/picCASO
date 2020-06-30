%The default results present are for E. coli model iAF1260

%% Gene information
load FBAresult.mat;
prompt1='Please input the serial number of the first gene from the FBAresult table that you would like to delete:\n';
serial1=input(prompt1);
prompt2='Please input the serial number of the second gene from the FBAresult table that you would like to delete:\n';
serial2=input(prompt2);
FBAresultCellArray=table2cell(FBAresult);
FBAresultArray=cell2mat(FBAresultCellArray(:,1));
ID1=find(FBAresultArray==serial1);
grRatio1=FBAresultCellArray{ID1,4};
ID2=find(FBAresultArray==serial2);
grRatio2=FBAresultCellArray{ID2,4};

%% Organism information
prompt="Please input the organism's generation time in minutes under the conditions being considered (on solid medium):\n";
GenTime=input(prompt);
ImprovedFitness = 10/GenTime; %Each step of the code translates to 10 minutes
BaseFitness1 = ImprovedFitness*grRatio1;
FitnessGain1 = ImprovedFitness - BaseFitness1;
BaseFitness2 = ImprovedFitness*grRatio2;
FitnessGain2 = ImprovedFitness - BaseFitness2;
save('Fitnesses.mat','BaseFitness1','FitnessGain1','BaseFitness2','FitnessGain2')

%% Device information
disp('Please pick AHL/device pairs using the software created by Kylilis et al., 2018.')
load('DeviceData.mat');
promptDevice1="Please choose the device for the first cell type out of the following.\nrhl | lux | tra | las | cin | rpa\n";
Device1=input(promptDevice1,'s');
promptAHL1="Please choose its AHL.\nC4 HSL | 3O C6 HSL | 30 C8 HSL | 3O C12 HSL | 3OH C14 | pC HSL\n";
AHL1=input(promptAHL1,'s');
promptDevice2="Please choose the device for the second cell type out of the following.\nrhl | lux | tra | las | cin | rpa\n";
Device2=input(promptDevice2,'s');
promptAHL2="Please choose its AHL.\nC4 HSL | 3O C6 HSL | 30 C8 HSL | 3O C12 HSL | 3OH C14 | pC HSL\n";
AHL2=input(promptAHL2,'s');
Match1 = find(string(DeviceData(:,1))==[Device1,' device']&string(DeviceData(:,2))==AHL1);
Threshold1=DeviceData{Match1,3};
Match2 = find(string(DeviceData(:,1))==[Device2,' device']&string(DeviceData(:,2))==AHL2);
Threshold2=DeviceData{Match2,3};
save('Thresholds.mat','Threshold1','Threshold2');