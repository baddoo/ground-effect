%   Returns a set of function handles corresponding to the circular arc wing
%   corresponding to annulus of interior radius q.
%   
%   Note that the length of the wing has not been normalised.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [f,fd,a,zt,d] = circularArcWing(gamma,q)

N = [1e3,1e2];
% Present the inital map
map = @(zVar) P(zVar.*conj(gamma),q,N)./(P(zVar./gamma,q,N).*P(conj(gamma),q,N)-P(1/gamma,q,N).*P(zVar.*conj(gamma),q,N));

% Differentiate initial map
mapd  = @(zVar) (Pd(zVar.*conj(gamma),q,N).*conj(gamma).*(P(zVar./gamma,q,N).*P(conj(gamma),q,N)-P(1/gamma,q,N).*P(zVar.*conj(gamma),q,N))- ...
                    P(zVar.*conj(gamma),q,N).*(Pd(zVar./gamma,q,N)./gamma.*P(conj(gamma),q,N)-P(1/gamma,q,N).*conj(gamma)*Pd(zVar.*conj(gamma),q,N)))...
                   ./(P(zVar./gamma,q,N).*P(conj(gamma),q,N)-P(1/gamma,q,N).*P(zVar.*conj(gamma),q,N)).^2;

% Find zeros of derivative on interior circle (correspond to pre-images locations of leading
% and trailing edges).
f = @(th) abs(mapd(q*exp(1i*th))); 
opts = optimset('Display','off');
thv(1) = fsolve(f,rand(1),opts);
% Now extract found root from function with smooth, periodic regulariser.
f1 = @(th) f(th)./abs(sin((th-thv(1))/2));
thv(2) = fsolve(f1,0,opts);        

% Sort zeros into leading and trailing edges
A1 = 1./exp(1i*angle(map(-1)-map(1i)));
[~,ind] = sort(real(A1*map(q.*exp(1i*thv))));
zt = q*exp(1i*thv(ind));

% Find length of plate to rescale
% zs = A1*map([zt,-q*exp(1i*mean(zt))]); 
% ma = imag(zs(2)-zs(1))./real(zs(2)-zs(1));
% mb = imag(zs(3)-zs(2))./real(zs(3)-zs(2));
% cx = (ma*mb*imag(zs(1)-zs(3)) + mb*real(zs(1)+zs(2)) - ma*real(zs(2)+zs(3)))...
%        ./(2*(mb-ma));
% cy = -1/ma*real(cx-(zs(1)+zs(2))/2) + imag(zs(1)+zs(2))/2;
% c = complex(cx,cy)
% rad = abs(zs(1) - c);
% (pi+angle(-(c-A1*map(zt(2)))./(c-A1*map(zt(1)))));

wingLength = 1;
A = A1./wingLength; s = -real(A*map(zt(1)))-1i*imag(A*map(-1));

% Find shift distance
f   = @(zVar) A*map(zVar) + s;
fd  = @(zVar) A*mapd(zVar);
%fdd = @(zVar) A*mapdd(zVar);

L = 1;
a = 1;

d = imag(f(zt(1)));

end