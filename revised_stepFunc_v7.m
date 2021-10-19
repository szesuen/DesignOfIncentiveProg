function [x, T, optClearanceRate] = revised_stepFunc_v7(params)

%original problem (step function)
% Version: JGv6_SS (SCRIPT), April 19, 2020
%
% min g(x)
% s.t.
% gamma*T <= B
% T(i) = -a(x(i)) - c(i)*x(i) +


PLOT_CHECK = 0;
% a = 780/17;  %780/17;
% c = 2*a*[5/26, 9/26, 13/26, 18/26, 22/26, 1]'; %c coefficient for each of the types
% xSQ_orig = c ./ 2 / a;
% N = length(c);  %number of types
% gamma = [0.0350242, 0.0434782, 0.0543478, 0.0772947, 0.0265701, 0.763285]';  %proportion in each type.  from data
% B = 0;

% Input parameters: rename to be ones from param
N = params.m;
c = params.c';
B = params.B;
gamma = params.gamma';
xSQ_orig = params.x_SQ;
a = params.a;

% Functions
aFunc =  @(x) (-a * pow_pos(x, 2));  %original quadratic version
aFuncPrime = @(x) (- 2 * a  * x);      %original quadratic version
nuFunc  = @(x) (gamma' *  (0.3891*log(  x  +0.0829) + 0.9690));  %version from original submission



% FIND STATUS-QUO LEVELS
xSQ = NaN(N, 1);
uSQ = NaN(N, 1);
for ii = 1:N
    xSQ(ii) = fzero(@(x) aFuncPrime(x) + c(ii), [0, 1]);
end
assert(sum(abs(xSQ - xSQ_orig)) <= 1E-7) %check, since we have x_SQ from input params



uSQ = aFunc(xSQ) + c.*xSQ;

phiFullFunc = @(j,x)  (uSQ(j)  - aFunc(x) - c(j)*x);


% check: 
if(PLOT_CHECK)
    t = linspace(0, 1, 1000);
    for ii = 1:N
        figure;
        plot(t, phiFullFunc(ii, t), 'b--');
    end
end
phiFunc = phiFullFunc;


sumVal = zeros(N,1);

cvx_begin
variable x(N,1)
variable T(N,1)
variable y(N, N)        % main trick is to use auxiliary variables to avoid the indicator.
expression sumVal(N, 1) % Use this to avoid double to cvx expression error
maximize( nuFunc(x) )
subject to

for ii = 1:N
    for jj = 1:ii
        if(jj == 1)  %c0 =0, x0=0.
            sumVal(ii) = phiFullFunc(jj, 0) ;
        else
            sumVal(ii) = sumVal(ii) + (c(jj) - c(jj-1))*x(jj-1) + phiFullFunc(jj, y(ii, jj));            
            y(ii, jj) <= x(jj - 1);
            y(ii, jj) <= xSQ(jj);
        end
    end

    % impose constraint
    T(ii) >= -aFunc(x(ii)) - c(ii)*x(ii) + sumVal(ii);
    T(ii) >= 0;  %Sze added, since T had neg values
end

%budget and utility constraints
% NOTE: ADDED Prime because assuming gamma is column vector (19 Apr 2020)
gamma' * T <= B;

%range constraints
xSQ <= x;
x <= ones(N,1) ;
x(2:N) >= x(1:N-1); % increasing
cvx_end


optClearanceRate = nuFunc(x);

% assert(sum(abs(x - xSQ_orig)) <= 1E-5) %check, since we have x_SQ from input params