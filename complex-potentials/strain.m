%   Returns function handles corresponding to the complex potential  and velocity
%   for a straining flow.
%   
%   strain(Q,N) evaluates the constant C for the annulus with interior radius Q 
%   and truncation N.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [potential,velocity] = strain(q,varargin)

switch nargin
    case 2
        N = [];
    case 3
        N = varargin{1};
end

coef = 1;

potential = @(zVar) coef*((zVar.*Pd(zVar,q,N)./P(zVar,q,N) ...
                   + zVar.^2.*Pdd(zVar,q,N)./P(zVar,q,N)...
                   - (zVar.*Pd(zVar,q,N)./P(zVar,q,N)).^2));

velocity = @(zVar)  coef*(Pd(zVar,q,N)./P(zVar,q,N) - 3*zVar.*(Pd(zVar,q,N)./P(zVar,q,N)).^2 ...
                    +2*zVar.^2.*(Pd(zVar,q,N)./P(zVar,q,N)).^3 + 3*zVar.*Pdd(zVar,q,N)./P(zVar,q,N)...
                    -3*zVar.^2.*Pd(zVar,q,N).*Pdd(zVar,q,N)./P(zVar,q,N).^2 + zVar.^2.*Pddd(zVar,q,N)./P(zVar,q,N));

end


