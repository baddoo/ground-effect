%   Returns function handles corresponding to the complex potential and complex velocity
%   for a uniform flow. Note that the complex velocity is the conjugate
%   of the actual velocity.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [potential,compVel] = uniform(q,a,varargin)

switch nargin
    case 2
        %N = [];
        N = [1e2,100];
    case 3
        N = varargin{1};
end

potential = @(zVar) a*zVar.*Pd(zVar,q,N)./P(zVar,q,N);

compVel = @(zVar)  a*(Pd(zVar,q,N)./P(zVar,q,N) +zVar.*Pdd(zVar,q,N)./P(zVar,q,N) ...
                    - zVar.*Pd(zVar,q,N).^2./P(zVar,q,N).^2);

end


