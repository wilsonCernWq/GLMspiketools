function [f,df,ddf] = logexp4(x)
%  [f,df,ddf] = logexp4(x);
%
%  Implements the nonlinearity:  
%    f(x)   = log(1+exp(x)).^pow  ==>  pow = 4
%    df(x)  = f'(x)
%    ddf(x) = f"(x)
%
%  General formulas:
%    f(x)  = log(1+exp(x))^k
%    f'(x) = k * log(1+e^x)^(k-1) * e^x/(1+e^x);
%    f"(x) = k(k-1) * log(1+e^x)^(k-2) * (e^x/(1+e^x))^2
%               + k * log(1+e^x)^(k-1) * e^x/(1+e^x)^2
%

pow = 4;
f0  = @(xx) logexp_pow(xx, pow);

if nargout <= 1
   f = f0(x);
elseif nargout == 2
   [f,df] = f0(x);
elseif nargout > 2
   [f,df,ddf] = f0(x);
end
