function [optimal_val_all, optimal_x_all, optimal_beta_all, optimal_That_all, optimal_T_all] = revised_test_convexFunc_JGv10(params)

% revised_test_convexFunc_JGv10
% -----------------------------------
% Extension of model (convex function)
% Version History 
% V4, April 30, 2020: basic version working for M = 0, ... N -1
% V5, May 6, 2020   : Version to calcualte all M simultaneously
% V6, May 12 2020   : debugging and testing
% V8, May 20 2020   : debugging
% V9, May 29 2020   : debugging - check for status "containing" Solved
%                     instead of exact string match.
%
% ASSUMPTION: a(x) = - a x^p  for p = 1 + 1/L, L = 1, 2, 3, .....
%


% For Debugging
DEBUG_CHECK = 0; % Set this to 1 to enable debugging

% Budget parameter
B = params.B;

xSQ = params.x_SQ; %[5/26, 9/26, 13/26, 18/26, 22/26, 1]'; % status quo utilities observed from data
gamma = params.gamma' ; %[0.0350242, 0.0434782, 0.0543478, 0.0772947, 0.0265701, 0.763285]';  %proportion in each type.  from data
N = params.m;  %number of types

% Utility Parameters
% NOTE: L = 1 is original quadratic case. For other values of L, we have to
% re-estimate a and c using the data. (not yet done)
L = 1;        
p = 1 + 1./L;   

% estimate a and c based on xSQ and utility function (see notes)
a = 60 ./ (p * (1 - xSQ(2)^(p-1)));
c = a * p * xSQ.^(p-1);

% check
if (L == 1)
    % case of our original manuscript
    a_orig = 780/17;  
    c_orig = 2*a*[5/26, 9/26, 13/26, 18/26, 22/26, 1]'; %c coefficient for each of the types
    if(norm(a - a_orig) + norm(c - c_orig) < 1E-7)
        disp('Check Passed: Estimates of a and c match that quadratic case');
    else 
        error('Check Passed: Estimates of a and c do not match quadratic case');
    end    
end


% --------------------------------------------------------------------------------
% Step 1: Generate polynomial coefficients for the beta * x(beta) - v(beta) terms
% NOTE: polynomial convention: increasing degree, i.e., start from 0, and then go to L
% -------------------------------------------------------------------------------
tmp = [1-p, 1];
for ii = 1:L
    tmp = conv([1, 1], tmp);
end
polybase1 = tmp;

% USAGE:
% testing = (c1.^(L+1:-1:0)) .* polybase1)
% constant factor = a /(ap)^(L+1)

% --------------------------------------------------------------------------------
% Step 2: Generate polynomial coefficients for the v_j(beta) - v_i(beta) terms
% --------------------------------------------------------------------------------
tmp = [1, 1];
for ii = 1:L
    tmp = conv([1, 1], tmp);
end
polybase2 = tmp;

% USAGE:
%c2 = c(2); c1 = c(1);
%testing = (c2.^(L+1:-1:0) - c1.^(L+1:-1:0)) .* polybase2.
% constant factor = a * (p - 1) /(ap)^(L+1)

% create constant to help simplify notation
POLYCONST = a ./ (a * p).^(L+1);

% Find Status-Quo Levels
uSQ_vec = POLYCONST * (p  - 1) * c.^(L+1);
uSQ1 = uSQ_vec(1);


% NOTE: We can't do M = N because we have a type that is fully adherent without any incentives. 
% Loop over M and repeatedly solve the model
optimal_val_all = NaN(1, N);
optimal_x_all   = NaN(N, N);
optimal_beta_all = NaN(N, N);
optimal_That_all = NaN(1, N);
optimal_T_all    = NaN(N, N);

% for debugging
optimal_Thatslack_all = NaN(1, N);
optimal_Tslack_all    = NaN(N, N);


% Range of M is from 0 to N-1
for M = 0:N-1
    % ---------------------------------------------------------------------
    % Begin CVX MODEL
    if(M == 0)
        % feasible case
        if(B >= uSQ1 + a - c(1))
            optimal_val  = 1;      % bug fix May 29, 2020
            optimal_x    = ones(N, 1);
            optimal_beta = zeros(N, 1);
            optimal_That= uSQ1 + a - c(1);
            optimal_T   = optimal_That * ones(N, 1);
            optimal_Thatslack = uSQ1 + a - c(1);
            optimal_Tslack    = optimal_Thatslack*ones(N, 1);            
        else
            % infeasible case
            optimal_val = -inf;
            optimal_x   = NaN(N, 1);
            optimal_beta = NaN(N, 1);
            optimal_That= NaN;
            optimal_T   = NaN(N, 1);
            optimal_Thatslack = NaN;
            optimal_Tslack    = NaN(N, 1);
        end        
    else
        % CASE where M > 0 
        cvx_begin
        variable mybeta(M,1)
        variable T(M,1)
        variable That
        expression vdiff_1step(M, 1)  % terms of this summand are v_{k+1} - v_k
        expression vdiff_Mstep(M, 1) % terms of this summand are v_{M+1} - v_i
        
        % Objective
        maximize(0.3242 * L * gamma(1:M)' * log((c(1:M) + mybeta)./(a * p)) + 1)
        subject to
        
        % build summand_vector to capture v_(k+1)  - v_k
        for ii = 1:M
            vdiff_1step(ii) = POLYCONST * (p - 1) * sum(((c(ii+1).^(L+1:-1:0) - c(ii).^(L+1:-1:0)) .* polybase2) .* [1, mybeta(ii), pow_pos(mybeta(ii), 2:(L+1))]);
        end
        
        % FIRST SET OF CONSTRAINTS
        for ii = 1:M
            % Added May 20, 2020 -- more explicit handling
            if(ii == 1)
                T(ii) >= uSQ1 + POLYCONST * sum((c(ii).^(L+1:-1:0)) .* polybase1 .* [1, mybeta(ii), pow_pos(mybeta(ii), 2:(L+1))]);
            else
                T(ii) >= uSQ1 + POLYCONST * sum((c(ii).^(L+1:-1:0)) .* polybase1 .* [1, mybeta(ii), pow_pos(mybeta(ii), 2:(L+1))]) + sum(vdiff_1step(1:ii-1));
            end
        end
        
        % SECOND SET OF CONSTRAINTS, only applies for M < N
        if(M < N)
            % build summand_vector to capture v_(M+1)  - v_i
            for ii = 1:M
                vdiff_Mstep(ii) = POLYCONST * (p - 1) * sum(((c(M+1).^(L+1:-1:0) - c(ii).^(L+1:-1:0)) .* polybase2) .* [1, mybeta(ii), pow_pos(mybeta(ii), 2:(L+1))]);
            end
            
            % Find Status-Quo Level for type M + 1
            uSQMplus1 = POLYCONST * (p  - 1) * c(M+1).^(L+1);
            
            % Bug fixed May 20, 2020
            % That >= uSQ1 + a - c(M+1) + uSQMplus1;    
            
            % Bug fixed Jun 29, 2020
            That >= uSQMplus1 + a - c(M+1);
            for ii = 1:M
                % Added May 20, 2020 -- more explicit handling
                if(ii == 1)
                    That >= uSQ1 + a - c(M+1) + vdiff_Mstep(ii);
                else
                    That >= uSQ1 + a - c(M+1) + vdiff_Mstep(ii) + sum(vdiff_1step(1:ii-1));
                end
            end
        end
        
        %budget constraint
        (1 - sum(gamma(1:M))) * That + gamma(1:M)' * T <= B;
        
        % T nonnegative 
        T >= 0;
        
        %range constraints
        mybeta >= 0; % added May 20, 2020 bugfix
        mybeta <= [mybeta(2:end); a * p  - c(M)];
        cvx_end
        
        
        % write out values
        %optimal_beta = [mybeta; ones(N - M, 1) * mybeta(end)];
        optimal_beta = [mybeta; NaN(N - M, 1)];
        optimal_Thatslack = That;
        optimal_Tslack    = [T; That* ones(N - M, 1)];
        
        % Remove slack from T
        optimal_T = NaN(M, 1);
        optimal_That = NaN;
        
        if(contains(cvx_status, 'Solved'))  % edited May 29, 2020 to "contains"
            for ii = 1:M
                optimal_T(ii) = uSQ1 + POLYCONST * sum((c(ii).^(L+1:-1:0)) .* polybase1 .* [1, mybeta(ii), pow_pos(mybeta(ii), 2:(L+1))]) + sum(vdiff_1step(1:ii-1));
            end
            optimal_That_vec = NaN(M, 1);
            for ii = 1:M
                optimal_That_vec(ii) = uSQ1 + a - c(M+1) + vdiff_Mstep(ii) + sum(vdiff_1step(1:ii-1));
            end
            optimal_That_vec = [0; optimal_That_vec];   % bugfix May 20, 2020
            optimal_That = max(optimal_That_vec);
            optimal_T = [optimal_T; optimal_That * ones(N - M, 1)];
        else
            % added May 20, 2020
            optimal_That = NaN;
            optimal_T = NaN(N, 1);
        end
        % END slack
        
        % edited May 29, 2020 to "contains"
        if(contains(cvx_status, 'Solved'))
            optimal_x = [((c(1:M) + mybeta)./(a * p)).^L; ones(N-M, 1)];
        else
            optimal_x = NaN(N, 1);
        end
        optimal_val = cvx_optval;       
    end
    
    % Write to output cells
    optimal_val_all(M+1) = optimal_val;
    optimal_x_all(:, M+1) = optimal_x;
    optimal_beta_all(:, M+1) = optimal_beta;
    optimal_T_all(:, M+1) = optimal_T;
    optimal_That_all(M+1) = optimal_That;
    
    % for debugging
    optimal_Tslack_all(:, M+1) = optimal_Tslack;
    optimal_Thatslack_all(M+1) = optimal_Thatslack;   

end

if(DEBUG_CHECK)
    % DEBUG And checking
    % check that beta and x are correctly related.
    check_mat_1 = ((repmat(c, 1, N) + optimal_beta_all) ./ (a * p)).^L - optimal_x_all;
    
    % construct Ti differently
    for M = 0:N-1
        optimal_beta = optimal_beta_all(:, M+1);
        optimal_x = optimal_x_all(:, M+1);
        vij = POLYCONST * (p-1) * (c + optimal_beta').^(L+1);
        test_T_vec = NaN(M, 1);
        for ii = 1:M
            test_T_vec(ii) = optimal_beta(ii) * optimal_x(ii);
            for kk = 1:ii
                if(kk == 1)
                    test_T_vec(ii) = test_T_vec(ii) + (uSQ1 - vij(kk, kk));
                else
                    test_T_vec(ii) = test_T_vec(ii) + (vij(kk, kk-1) - vij(kk, kk));
                end
            end
        end
        diff_T = optimal_T_all(1:M, M+1) - test_T_vec;
    end
    
    % construct 1 step v differences
    dv = zeros(M, 1);
    cum_dv = zeros(M, 1);
    for ii = 1:M
        cum_dv(ii) = sum(dv);
        dv(ii) = vij(ii+1, ii) - vij(ii, ii);
    end
end