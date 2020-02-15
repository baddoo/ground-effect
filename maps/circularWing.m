%   Returns a set of function handles corresponding to the circular wing
%   map of unit length corresponding to annulus of interior radius q.
%   
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [f,fd,a,d] = circularWing(q)

f   = @(zVar) (1-q.^2)./(2i*q).*(zVar+1)./(zVar-1);
fd  = @(zVar) (1-q.^2)./(2i*q).*(-2)./(zVar-1).^2;
%fdd = @(zVar) (1-q.^2)./(2i*q).*(4)./(zVar-1).^3;

a = (1-q.^2)./(1i*q);
d = (q+q^-1)/2;
end