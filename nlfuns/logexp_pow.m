function [f,df,ddf] = logexp_pow(x,pow)
%  [f,df,ddf] = logexp_pow(x);
%
%  Implements the nonlinearity:  
%     f(x) = log(1+exp(x)).^pow
%  plus first and second derivatives
%
%  General formulas:
%    f(x)  = log(1+exp(x))^k
%    f'(x) = k * log(1+e^x)^(k-1) * e^x/(1+e^x);
%    f"(x) = k(k-1) * log(1+e^x)^(k-2) * (e^x/(1+e^x))^2
%               + k * log(1+e^x)^(k-1) * e^x/(1+e^x)^2
%


f0 = log(1+exp(x));
f = f0.^pow;

if nargout > 1
    df = pow*f0.^(pow-1).*exp(x)./(1+exp(x));
end
if nargout > 2
    if pow == 1
        ddf = pow*f0.^(pow-1).*exp(x)./(1+exp(x)).^2;
    else
        ddf = pow*f0.^(pow-1).*exp(x)./(1+exp(x)).^2 + ...
              pow*(pow-1)*f0.^(pow-2).*(exp(x)./(1+exp(x))).^2;
    end
end
