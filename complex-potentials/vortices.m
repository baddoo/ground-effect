%   Returns function handles corresponding to the complex potential and  complex velocity
%   for a number of vortices inside the annulus. Note that the complex velocity is the conjugate
%   of the actual velocity.
%   
%   vortices(BETA,Q,N) evaluates the constant C for the annulus with interior radius Q 
%   and truncation N.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [potential,compVel] = vortices(beta,q,varargin)

switch nargin
    case 2
        N = [];
    case 3
        N = varargin{1};
end
    
beta3 = permute(beta(:),[3,2,1]);

potential = @(zVar) 1/(2i*pi)*sum(log(abs(beta3)*P(zVar./beta3,q,varargin{1})...
                                        ./P(zVar*conj(beta3),q,varargin{1})),3);
compVel = @(zVar) 1/(2i*pi)*sum(1./beta3*Pd(zVar./beta3,q,N)./P(zVar./beta3,q,N)...
                    - conj(beta3)*Pd(zVar*conj(beta3),q,N)./P(zVar.*conj(beta3),q,N),3);

end