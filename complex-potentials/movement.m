%   Returns function handles corresponding to the complex potential  and velocity
%   for the movement of the wing.
%   
%   vortices(BETA,Q,N) evaluates the constant C for the annulus with interior radius Q 
%   and truncation N.
%
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function [potential,compVel] = movement(f,d,dDot,aDot,q,varargin)

switch nargin
    case 5
        N = [];
    case 6
        N = varargin{1};
end

% Number of points in integration. Should vary depending on application.
ni = 200;

tInt = linspace(0,2*pi,ni + 1); tInt(end) = [];
tInt4 = permute(tInt,[1,4,3,2]);
zInt4 = q*exp(1i*tInt4); 

Mprime = @(zVar) imag(-1i*dDot*(f(zVar)+1i*d) + 1i*aDot/2*abs(f(zVar)+1i*d).^2);

% Integrate using exponentially convergent trapeziod rule

I = -1/ni*sum(Mprime(zInt4),4);

M = @(zVar) Mprime(zVar) + I;

% Rearrange to avoid errors close to interior radius
g = @(zVar) M(zVar)./zVar.^2;

% Regularise kernel of potential function
regKerPot = @(zVar0,zVar1) Pd(zVar0./zVar1,q,N)./P(zVar0./zVar1,q,N) + zVar1./(zVar1-zVar0);
pot1 = @(zVar)  zVar/pi.*sum(1i*zInt4.*g(zInt4).*regKerPot(zVar,zInt4),4)*2*pi/ni;
pot2 = @(zVar) -zVar/pi.*sum(1i*zInt4.^2.*(g(zInt4)-g(zVar))./(zInt4-zVar),4)*2*pi/ni;

% Combine potentials
potential = @(zVar) pot1(zVar) + pot2(zVar);

% Now compute velocity (can regularise in future)

% Define kernel of 
velKer =  @(zVar0,zVar1) Pd(zVar0./zVar1,q,N)./P(zVar0./zVar1,q,N) ...
         + zVar0./zVar1.*Pdd(zVar0./zVar1,q,N)./P(zVar0./zVar1,q,N) ...
         - zVar0./zVar1.*Pd(zVar0./zVar1,q,N).^2./P(zVar0./zVar1,q,N).^2;
compVel = @(zVar) 1/pi.*sum(1i*zInt4.*g(zInt4).*velKer(zVar,zInt4),4)*2*pi/ni;

end