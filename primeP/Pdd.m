%   Evaluates the second derivative of the prime function P.
%   
%   Pd(ZETA,Q,N) evaluates the derivative of the prime function Pd at the 
%   points ZETA for the annulus with interior radius Q and truncation 
%   vector N. The first entry of N corresponds to the truncation for the 
%   constant C and the second entry of N corresponds to the truncation of 
%   the rest of the series.
%   
%   ZETA can be a matrix, but Q must be a scalar.
%
%   Q is 1/2 if left empty.
%
%   If left empty, N is selected to be such that the truncation error is below a
%   tolerance of 1e-3, and 1e3 terms are used to compute C.
%
%   This function uses the Laurent expansion for the prime function P. The
%   truncation error is bounded by a prescribed tolerance using the
%   integral test.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function Pdd = Pdd(zeta,varargin)

% This tells us the error for truncation nS
err =  @(nS,qv,nC) abs(C(qv,nC)*qv/2/log(1/qv).^1.5)...
                *(sqrt(pi)*erfc((nS-1)*sqrt(log(1/qv)))...
                  + 2*nS*qv^((nS-1)^2)*sqrt(log(1/qv)));
        

% Selection of N for non-derivative to be used as initial guess
bound = @(vtol,qv,nC) round(1 + erfcinv(qv*vtol/abs(C(qv,nC))*sqrt(log(1/qv)/pi))./sqrt(log(1/qv)));

opts = optimset('Display','off');

switch nargin
    case 1
        q = 1/2;
        tol = 1e-3;
        N(1) = 1e3;
        N(2) = round(fsolve(@(nS) tol - err(nS,q,N(1)),bound(tol,q,N(1)),opts));
    case 2
        q = varargin{1};
        tol = 1e-3;
        N(1) = 1e3;
        N(2) = round(fsolve(@(nS) tol - err(nS,q,N(1)),bound(tol,q,N(1)),opts));
    case 3
        q = varargin{1};
        N = varargin{2};
        if isempty(N)
            tol = 1e-3;
            N(1) = 1e3;
            N(2) = round(fsolve(@(nS) tol - err(nS,q,N(1)),bound(tol,q,N(1)),opts));
        end  
end

n2 = permute(-N(2):N(2),[1,3,2]);
Pdd = C(q,N(1))*sum((-1).^n2.*q.^(n2.*(n2-1)).*n2.*(n2-1).*zeta.^(n2-2),3);

end