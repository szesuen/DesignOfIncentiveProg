clear
clear all
close all
clc


%% Get inputs from input file
timeStamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');
parentDir = cd;
inputs2;


%% Different Optimization Cases
for incType = 1:3 %1=linear, 2=step, 3=convex
    outMatAll{incType} = [];  %Array to hold output
    i = 1; 
    
    for WVal = 0:0.5:35  %Sweep over different budget levels
        params.B = WVal;
        params.useStepFunc = incType-1;  %0 or 1 or 2.  0: linear incentives, 1: step func, 2: convex
        
        % Sweep over different utility functional forms
        for nonConvUtilFuncIdx = 1 %1:5 %1=only quadratic case.  1:5  %only changes anything if find_optX_v6.m is using the non-convex solver
            params.utilityFuncType = nonConvUtilFuncIdx;  %1=Polynomial, 2=Logarithmic, 3=Exponential, 4=Negative Powers, 5=Fractional Powers
            
            %optimize
            find_optX_v6
            
            %save output
            outMatAll{incType}(i,:) = optX';
            i = i+1;
        end
    end
end

figureFolderNameStr = strcat('figures','_bud',num2str(params.budgetVersion),'_', timeStamp);
mkdir(figureFolderNameStr)
cd(figureFolderNameStr)
save outputsFile
cd(parentDir)
