%   Returns function handles corresponding to the complex potential and complex velocity
%   for the (clockwise) circulation around the wing.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [potential,compVel] = circulatory

potential = @(zVar) -1/(2i*pi)*log(zVar);                                        
compVel = @(zVar) -1/(2i*pi*zVar);

end