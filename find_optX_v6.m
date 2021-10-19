%% optimize using different incentive strucures
if (params.useStepFunc == 0)  %linear incentives for comparison
    optX = min((params.B+params.c)'/(2*params.a),1);
    
elseif (params.useStepFunc == 1)  %% step function
    [optX, optT, optClearRate] = revised_stepFunc_v7(params);

elseif (params.useStepFunc == 2)  %% Convex incentives (piecewise linear (PW) incentives)
    
    %%%% SPECIFY INTEGER L IN revised_test_convexFunc_JGv10.m 
    [optimal_val_all, optimal_x_all, optimal_beta_all, optimal_That_all, optimal_T_all] = ...
        revised_test_convexFunc_v10(params);
    
    %gives values across all M, need to choose the correct value
    notValidIdx = (optimal_val_all == -Inf);
    optimal_val_all(notValidIdx) = NaN;
    [optVal, optIdx] = max(optimal_val_all);
    
    optX = optimal_x_all(:,optIdx);
    optBeta = optimal_beta_all(:,optIdx);
    optT_hat = optimal_That_all(:,optIdx);
    optT = optimal_T_all(:,optIdx);

else
    disp('params.useStepFunc is not 0 or 1 or 2');
end

