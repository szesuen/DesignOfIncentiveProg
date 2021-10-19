function dxdt = dynamics(t,x,params)
% v2struct(params);

switch params.optCase
  case 0
    nu = params.nuSQ;
  case 1
    nu = params.nuOpt;
end

S = x(1);
E = x(2);
I = x(3);
N = sum(x);

dxdt(1,1) = -params.beta*(I/N)*S + (nu*I) - (params.d*S)+ params.b*N;
dxdt(2,1) = params.beta*(I/N)*S - params.alpha*E - params.d*E;
dxdt(3,1) = params.alpha*E-nu*I - params.mu*I;

end