%   Evaluates the prime function P.
%   
%   P(ZETA,Q,N) evaluates the prime function P at the points ZETA for the
%   annulus with interior radius Q and truncation vector N. The first entry
%   of N corresponds to the truncation for the constant C and the second
%   entry of N corresponds to the truncation of the rest of the series.
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


function P = P(zeta,varargin)

% This function tells you what N you should use to truncate the series to
% achieve a tolerance of vtol for interior radius qv and C truncation nC.
bound = @(vtol,qv,nC) round(1 + erfcinv(qv*vtol/abs(C(qv,nC))*sqrt(log(1/qv)/pi))./sqrt(log(1/qv)));

switch nargin
    case 1
        q = 1/2;
        tol = 1e-3;
        N(1) = 1e3;
        N(2) = bound(tol,q,N(1));
    case 2
        q = varargin{1};
        tol = 1e-3;
        N(1) = 1e3;
        N(2) = bound(tol,q,N(1));
    case 3
        q = varargin{1};
        N = varargin{2};
        if isempty(N)
            tol = 1e-3;
            N(1) = 1e3;
            N(2) = bound(tol,q,N(1));  
        end  
end

n2 = permute(-N(2):N(2),[1,3,2]);
P = C(q,N(1))*sum((-1).^n2.*q.^(n2.*(n2-1)).*zeta.^n2,3);

end