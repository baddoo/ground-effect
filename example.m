addpath('maps','primeP','complex-potentials')

q = .2;
nr = 100; nt = 100;
zr = linspace(q,1,nr+2); %zr(1) = []; zr(end) = [];
zt = linspace(0,2*pi,nt);
zeta = zr.*exp(1i*zt');
zetab = [q,1].*exp(1i*zt');
plot(zeta);

[f,fd,a] = circularWing(q);

flow = 'strain';
%flow = 'vortices';
flow = 'uniform';

if strcmp(flow,'strain')

[potential,compVel] = strain(q,a);
vals = imag(potential(exp(1i*[+1;-1]*pi/4).*linspace(q,1,10)));
clim = [0,3];

elseif strcmp(flow,'vortices')
    
    beta = -sqrt(q);
[potential,compVel] = vortices(beta,q,[1e3,100]);
vals = imag(potential(exp(1i*[+1;-1]*pi/4).*linspace(q,1,10)));
clim = [0,1];

elseif strcmp(flow,'uniform')
    
[potential,compVel] = uniform(q,a);
vals = imag(linspace(potential(-1),potential(.9),100));

clim = [0,1.5];

end

z = f(zeta);
zb = f(zetab);

ax = [-4,4,0,15];
subplot(2,1,1)
plot(z,'k','LineWidth',2);
hold on
plot(zb(:,1),'r','LineWidth',3)
plot(zb(:,2),'b','LineWidth',3)
hold off
axis equal
axis(ax)

subplot(2,1,2)
pcolor(real(z),imag(z),real(compVel(zeta)./fd(zeta)));
colormap jet; shading interp; caxis(clim); colorbar
hold on
contour(real(z),imag(z),imag(potential(zeta)),vals(:),'w','LineWidth',1);
plot(zb(:,1),'r','LineWidth',3)
plot(zb(:,2),'b','LineWidth',3)
hold off
axis equal
axis(ax)
return
%% 

disp('Displays flow.')
pause

%%

disp('Displays another flow');
pause