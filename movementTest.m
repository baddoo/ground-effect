
nInt = 1e2;
th = linspace(0,2*pi,nInt+1); th(end) = [];

f = @(zVar) real(zVar);%./zVar.^2;

q = .5;
zeta = sqrt(q);
integrand = @(zVar) f(zVar).*Pd(zeta./zVar,q)./P(zeta./zVar,q);
exact = integral(@(thv) 1i*q*exp(1i*thv).*integrand(q*exp(1i*thv)),0,2*pi)

zInt = q*exp(1i*th);
k3 = permute(1:1e2,[1,3,2]);
zt = linspace(q,1);
pI1 = Pd(zt,q)./P(zt,q);

pI2 =  1./(zt-1)...
          + sum(1./(zt-q.^(-2*k3))  - 1./zt + 1./(zt-q.^(2*k3)),3);
     
%plot(imag(pI1)); hold on; plot(imag(pI2),'--'); hold off;
      
partialIntegrand = zInt./(zeta-zInt)...
          + sum(zeta./(zInt-zeta.*q.^(-2*k3)) - zeta./zInt + zeta./(zInt - zeta.*q.^(2*k3)),3);
          
%partialIntegrand = -1 - zeta./(zInt-zeta)...
%          + sum(-q.^(2*k3) - q.^(-2*k3) - zInt./zeta...
%          - zeta*(q.^(4*k3)./(zInt-zeta.*q.^(-2*k3)) + q.^(-4*k3)./(zInt-zeta.*q.^(2*k3))),3);
partialIntegrand = -1 - zeta./(zInt-zeta)...
          + sum(-q.^(2*k3)+(zeta.*q.^(4*k3))./(zeta.*q.^(2*k3)-zInt) - zInt./(zeta) + (zInt)./(zeta - zInt.*q.^(2*k3)),3);    

%plot(real(partialIntegrand)); hold on; plot(real(Pd(zeta./zInt,q)./P(zeta./zInt,q)),'--'); hold off;
regKer = @(ziv,zv) Pd(zv./ziv,q,[1e3,1e3])./P(zv./ziv,q,[1e3,1e3])-ziv./(zv-ziv);% zv./(ziv-zv);
trap = sum(1i*zInt.*integrand(zInt))*2*pi/nInt
trap = sum(1i*zInt.*f(zInt).*regKer(zInt,zeta))*2*pi/nInt% + zeta.*(f(zInt))./(zInt-zeta))*2*pi/nInt