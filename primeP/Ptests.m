N = 50;
Ninf = 100;
nq = 5;

q = linspace(.9,1,nq+2); q(1) = []; q(end) = [];
circ=exp(1i*linspace(0,2*pi,200));

err = zeros(N,nq);
bound = zeros(N,nq);
errP = zeros(1,nq);
for l = 1:nq
    
    pts = [q(l)*circ,1/q(l)*circ];
    
for j = 1:N
    
err(j,l) = norm(PfunDD(pts,q(l),j')-PfunDD(pts,q(l),Ninf),'inf');
%err(j,l) = norm(PfunD(pts,q(l),j')-PfunD(pts,q(l),Ninf),'inf');

PHI = @(nv,zv) exp(-.25*log(abs(zv)).^2./log(q(l)))...
             .*(1+erf(((2*nv-1)*log(q(l))+log(abs(zv)))./(2*sqrt(-log(q(l))))))...
             .*sqrt(-pi*abs(zv)./log(q(l)))/2./q(l)^.25;
% Zeroth derivative         
PHI = @(nv,zv) sqrt(pi/log(1/q(l)))*erfc((nv-1)*sqrt(log(1/q(l))))/q(l); 
% First derivative
PHI = @(nv) 1e-3-(q(l)^((nv-1)^2)+sqrt(pi*log(1/q(l)))*erfc((nv-1)*sqrt(log(1/q(l)))))...
                /q(l)/log(q(l));
% Second derivative            
PHI = @(nv) abs(C(q(l),1e3)*q(l)/2/log(1/q(l)).^1.5)...
                *(sqrt(pi)*erfc((nv-1)*sqrt(log(1/q(l))))...
                  + 2*nv*q(l)^((nv-1)^2)*sqrt(log(1/q(l))));
        
%tic;a = fsolve(PHI,100); toc;
bound(j,l) = abs(calcA(q(l),j))*norm(PHI(j),'inf');         

end

%errP(l) = norm(PfunAd(pts,q(l))-Pfun(pts,q(l),Ninf),'inf');
%errP(l) = norm(PfunDAd(pts,q(l))-Pfun(pts,q(l),Ninf),'inf');


end

%% Plots

cols = hot(floor(nq*2));
figure(1); clf;
for l = 1:nq
semilogy(err(:,l),'--o','Color',cols(l,:),'LineWidth',3);
hold on
semilogy(bound(:,l),'-x','Color',cols(l,:),'LineWidth',3);
end
ylim([1e-16,1e2])
grid on
hold off

function A = calcA(q,N)

n1 = permute(1:1e3,[1,3,2]);
An = prod((1+q.^(2*n1)).^2,3);
Ad = sum(q.^(n1.*(n1-1)),3);
A = An./Ad;

end

function P = PfunAd(zVar,q)

tol = 1e-3;
N = round(1 + erfcinv(q*tol/abs(calcA(q))*sqrt(log(1/q)/pi))./sqrt(log(1/q)));
              
n2 = permute(-N:N,[1,3,2]);
n1 = permute(1:1e3,[1,3,2]);
An = prod((1+q.^(2*n1)).^2,3);
Ad = sum(q.^(n1.*(n1-1)),3);
A = An./Ad;

P = A*sum((-1).^n2.*q.^(n2.*(n2-1)).*zVar.^n2,3);

end


function P = Pfun(zVar,q,N)

n2 = permute(-N:N,[1,3,2]);
n1 = permute(1:1e3,[1,3,2]);
An = prod((1+q.^(2*n1)).^2,3);
Ad = sum(q.^(n1.*(n1-1)),3);
A = An./Ad;

P = A*sum((-1).^n2.*q.^(n2.*(n2-1)).*zVar.^n2,3);

end

function Pp = PfunProd(zVar,q,N)

n1 = permute( 1:N,[1,3,2]);
Pp = (1-zVar).*prod((1-q.^(2*n1).*zVar).*(1-q.^(2*n1)./zVar),3);

end

function Pd = PfunD(zVar,q,N)

n2 = permute(-N:N,[1,3,2]);
n1 = permute( 1:1e3,[1,3,2]);
An = prod((1+q.^(2*n1)).^2,3);
Ad = sum(q.^(n1.*(n1-1)),3);
A = An./Ad;

Pd = A*sum((-1).^n2.*q.^(n2.*(n2-1)).*n2.*zVar.^(n2-1),3);

end

function Pdd = PfunDD(zVar,q,N)

n2 = permute(-N:N,[1,3,2]);
n1 = permute( 1:1e3,[1,3,2]);
An = prod((1+q.^(2*n1)).^2,3);
Ad = sum(q.^(n1.*(n1-1)),3);
A = An./Ad;

Pdd = A*sum((-1).^n2.*q.^(n2.*(n2-1)).*n2.*(n2-1).*zVar.^(n2-2),3);

end