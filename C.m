%   Evaluates the constant C.
%   
%   C(Q,N) evaluates the constant C for the annulus with interior radius Q 
%   and truncation N.
%
%   If left empty, N is selected be 1e3.
%
%   Cite: Exact solutions for ground effect, P. J. Baddoo, M. Kurt, L. J.
%         Ayton, K. W. Moored, JFM Rapids, 2020


function C = C(q,N)

n1 = 1:N;
Cn = prod((1+q.^(2*n1)).^2,2); % Numerator of C
Cd = sum(q.^(n1.*(n1-1)),2); % Denominator of C
C = Cn./Cd;

end