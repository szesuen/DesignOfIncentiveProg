%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% inputs %%%%%%%%
 
%optimization problem vars
r     = 0.03; %0.5; 0.03; %discount factor
params.objectivetype = 2;  %log objective (use 1 for linear)  
params.budgetVersion = 0; %if nonzero, uses the budget constrained version of the problem (instead of NMB objective)
W     = 3*1450*(0.66); %willingness to pay factor

%utility function parameters
QALYvec = [1 1 0.6];

% a = 7.65;  %utility func sq term coefficient
% a = ((60*(9/26))/(1-(9/26)))  * (26/9) * (1/2);  %assuming linear reward

%using 60 = T' (see mathematica file)
a = 780/17;  %780/17;

c = 2*a*[5/26, 9/26, 13/26, 18/26, 22/26, 1]; %c coefficient for each of the types
m = length(c);  %number of types
x_SQ = c'/(2*a); %[0.192, 0.346, 0.500, 0.692, 0.846, 1 ]';

gamma = [0.0350242, 0.0434782, 0.0543478, 0.0772947, 0.0265701, 0.763285];  %proportion in each type.  from data

% epi vars
N_0 = 100000;
S_0 = N_0*(0.6-0.00195); %1000; % N_0*0.6; %10000;
E_0 = N_0*0.44; %50; %N_0*0.3975; %1000;  %S_0*(0.6/0.4);
I_0 = N_0*0.00195; %N_0*0.0021; %1000 ; %N_0*0.0025;  %0.0025*(S_0+E_0);  %0.0025*(S_0+E_0);
% N_0 = S_0 + E_0 + I_0;



beta  = 12.084866; %12.084589; %12.084866 ; %5.15; 4.65; %3.04; %2;% 10; %2; %transmission rate (10-15 ppl in 3 yrs)
alpha = 0.00479545454; %0.005275; % N*0.4*params.alpha = 211*N/100000, where 211 from inci rate WHO ....   or alternatively %Trauer, 5-15% lifetime activation rate (67.4 yrs expectation of life) = 0.017 to 0.054 using p = 1-exp(-r*t) -- r = -ln(1-p)/t 

d     = 0.006718;   %background mortality rate (population average using India 2011 Census and WHO lifetables)
b     = d+0.011; % 0.00771; %birth rate 
mu    = -log(1-(420/2581.8)) + d ;    %TB mortality rate. 
trtCovRate = -log(1-(1936.158/2581.8));  % Number of notified cases 
nuSQ = trtCovRate * (gamma*(0.3891*log(  min( x_SQ,1)  +0.0829)+0.9690)) ; 

% %%coefficients for objective function
% For Convex Version
coef = zeros(m,m);
for i = 1:m
    for j = 1:i
        coef(i,j) = 2*gamma(i)*(c(i)-c(j));
    end
end
QcoeffMat = zeros(m,m);
for i = 1:m
    for j = 1:m
        QcoeffMat(i,j) = sum(gamma(max(i,j):m));
    end
end

% For Step Function Version
cDiff = c - [0,c(1:end-1)];
cDiffMat = tril(ones(m,m)).*repmat(cDiff,m,1);


integralCoef = trtCovRate/r ;




%% Gather all input parameters into a structure
optCase = 0;  %switch to one if you want nu to be nuOpt
params.r = r;
params.W = W;
params.a = a;
params.c = c;
params.gamma = gamma;
params.m = m;
params.x_SQ = x_SQ;
params.coef = coef;
params.QcoeffMat = QcoeffMat;
params.cDiffMat = cDiffMat;
params.N_0 = N_0;
params.S_0 = S_0;
params.E_0 = E_0;
params.I_0 = I_0;
params.beta = beta;
params.alpha = alpha;
params.b = b;
params.d = d;
params.mu = mu;
params.trtCovRate = trtCovRate;
params.nuSQ = nuSQ;
params.optCase = optCase;


%% Check calibration
initial_x = [S_0; E_0; I_0];
tspan = (0:1:100);
params.optCase = 0;
[t,xSQ] = ode45(@(t,y) dynamics(t,y,params), tspan, initial_x);

tol = 0.001;
prev = xSQ./repmat(sum(xSQ,2),1,3);
if (abs(prev(1,1) - prev(end,1)) > tol)    
    beta = betaCalibration(params, tol);
    params.beta = beta;
    fprintf('Beta recalibrated to %f \n', beta);
end

