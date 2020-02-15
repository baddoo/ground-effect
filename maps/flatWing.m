%   Returns a set of function handles corresponding to the flat plate wing
%   map of unit length corresponding to annulus of interior radius q.
%   
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [f,fd,a,zt,d] = flatWing(alpha,q)

N = [1e3,1e2];
% Present the inital map
if alpha ~= 0 
map = @(zVar) P(zVar.*exp(2i*alpha),q,N)./P(zVar,q,N);

% Differentiate initial map
mapd  = @(zVar) exp(2i*alpha)*Pd(zVar.*exp(2i*alpha),q,N)./P(zVar,q,N) ...
                  -P(zVar.*exp(2i*alpha),q,N).*Pd(zVar,q,N)./P(zVar,q,N).^2;
%mapdd = @(zVar)   -2*exp(2i*alpha)*Pd(zVar.*exp(2i*alpha),q,N).*Pd(zVar,q,N)./P(zVar,q,N).^2 ...
%                  + exp(4i*alpha)*Pdd(zVar,q,N)./P(zVar,q,N)...
%                  +P(zVar.*exp(2i*alpha),q,N).*(2*Pd(zVar,q,N).^2./P(zVar,q,N).^3 - Pdd(zVar,q,N)./P(zVar,q,N).^2);

% Find zeros of derivative on interior circle (correspond to pre-images locations of leading
% and trailing edges).
f = @(th) abs(mapd(q*exp(1i*th))); 
%fp = @(th) real(mapd(q*exp(1i*th)).*conj(1i*q*exp(1i*th).*mapdd(q*exp(1i*th))))./abs(mapd(q*exp(1i*th)));
% Set tolerance for root finding
% tol = 1e-10;
opts = optimset('Display','off');
thv(1) = fsolve(f,0,opts);
%thv(1) = mod(myNewton(f,fp,0,tol),2*pi);
% Now extract found root from function with smooth, periodic regulariser.
f1 = @(th) f(th)./abs(sin((th-thv(1))/2));
% f1p= @(th) (fp(th).*abs(sin((th-thv(1))/2))-f(th).*sign(sin((th-thv(1))/2))/2.*cos((th-thv(1))/2))...
%             ./abs(sin((th-thv(1))/2)).^2;
thv(2) = fsolve(f1,0,opts);        
%thv(2) = myNewton(f1,f1p,0,tol);

% Sort zeros into leading and trailing edges
A1 = sign(-alpha)*exp(-1i*alpha);
[~,ind] = sort(real(A1*map(q.*exp(1i*thv))));
zt = q*exp(1i*thv(ind));

% Find length of plate to rescale
plateLength = abs(diff(A1*map(zt)));
A =  A1./plateLength;

k3 = permute(1:1e3,[1,3,2]);
L = prod((1-q.^(2*k3)).^2,3);
a = -A*P(exp(2i*alpha),q,N)./L;

else

    
A1 = -1i;    
map = @(zVar) zVar.*Pd(zVar,q,N)./P(zVar,q,N);

% Differentiate initial map
mapd  = @(zVar) Pd(zVar,q,N)./P(zVar,q,N) + zVar.*Pdd(zVar,q,N)./P(zVar,q,N)...
                -zVar.*Pd(zVar,q,N).^2./P(zVar,q,N).^2;
% Find zeros of derivative on interior circle (correspond to pre-images locations of leading
% and trailing edges).
f = @(th) abs(mapd(q*exp(1i*th)));

opts = optimset('Display','off');
thv(1) = fsolve(f,pi/4,opts);
f1 = @(th) f(th)./abs(sin((th-thv(1))/2));
thv(2) = fsolve(f1,pi/4,opts);        

% Sort zeros into leading and trailing edges
[~,ind] = sort(real(A1*map(q.*exp(1i*thv))));
zt = q*exp(1i*thv(ind));

plateLength = abs(diff(A1*map(zt)));
A =  A1./plateLength;

a = A;

end

s = -real(A*map(zt(1))) - 1i*imag(A*map(-1));
f   = @(zVar) A*map(zVar) + s;
fd  = @(zVar) A*mapd(zVar);
d = imag(f(zt(1)));    

