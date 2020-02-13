%   Returns a set of function handles corresponding to the circular arc wing
%   map of unit length corresponding to annulus of interior radius q.
%   
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [f,fd,a,zt,d] = centeredCircularArcWing(phi,q)

N = [];
% Present the inital map
map = @(zVar) -P(sqrt(zVar).*exp(1i*phi),q,N).*P(-sqrt(zVar).*exp(1i*phi),q,N).*P(q*sqrt(zVar),q,N).*P(-q*sqrt(zVar),q,N)./...
              (P(q*sqrt(zVar).*exp(1i*phi),q,N).*P(-q*sqrt(zVar).*exp(1i*phi),q,N).*P(sqrt(zVar),q,N).*P(-sqrt(zVar),q,N)).*exp(-1i*phi);

% Differentiate initial map
mapd  = @(zVar) -((Pd(zVar*exp(2i*phi),q.^2,N)*exp(2i*phi).*P(q^2*zVar,q^2,N)+P(zVar*exp(2i*phi),q.^2,N).*Pd(q^2*zVar,q^2,N)*q^2).*...
                  P(q^2*zVar.*exp(2i*phi),q^2,N).*P(zVar,q^2,N)...
                  -(q^2*exp(2i*phi)*Pd(q^2*zVar.*exp(2i*phi),q^2,N).*P(zVar,q^2,N)+P(q^2*zVar.*exp(2i*phi),q^2,N).*Pd(zVar,q^2,N)).*...
                  P(zVar*exp(2i*phi),q.^2,N).*P(q^2*zVar,q^2,N))...
                  ./(P(q^2*zVar.*exp(2i*phi),q^2,N).*P(zVar,q^2,N)).^2.*exp(-1i*phi);        

% Find zeros of derivative on interior circle (correspond to pre-images locations of leading
% and trailing edges).
f = @(th) abs(mapd(q*exp(1i*th))); 
%fp = @(th) real(mapd(q*exp(1i*th)).*conj(1i*q*exp(1i*th).*mapdd(q*exp(1i*th))))./abs(mapd(q*exp(1i*th)));
% Set tolerance for root finding
% tol = 1e-10;
opts = optimset('Display','off');
thv(1) = fsolve(f,rand(1),opts);
%thv(1) = mod(myNewton(f,fp,0,tol),2*pi);
% Now extract found root from function with smooth, periodic regulariser.
f1 = @(th) f(th)./abs(sin((th-thv(1))/2));
% f1p= @(th) (fp(th).*abs(sin((th-thv(1))/2))-f(th).*sign(sin((th-thv(1))/2))/2.*cos((th-thv(1))/2))...
%             ./abs(sin((th-thv(1))/2)).^2;
thv(2) = fsolve(f1,0,opts);        
%thv(2) = myNewton(f1,f1p,0,tol);

% Sort zeros into leading and trailing edges
A1 = 1./exp(1i*angle(map(-1)-map(1i)));
[~,ind] = sort(real(A1*map(q.*exp(1i*thv))));
zt = q*exp(1i*thv(ind));

% Find length of plate to rescale
wingLength = angle((A1*map(zt(1)))./(A1*map(zt(2))));
A = A1./wingLength; s = -real(A*map(zt(1)))-1i*imag(A*map(-1));

% Find shift distance
f   = @(zVar) A*map(zVar) + s;
fd  = @(zVar) A*mapd(zVar);
%fdd = @(zVar) A*mapdd(zVar);

L = 1;
a = 1;

d = imag(f(zt(1)));

end