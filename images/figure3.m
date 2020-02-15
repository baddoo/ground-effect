addpath('../maps','../primeP','../complex-potentials')
q = 0.3;
%%%%%%%%%%%%%% Define mappings

beta = 1.5*q*exp(1i*(pi/4+.2));
phi = 5*pi/8;
c0 = exp(1i*linspace(0,2*pi)');
c1 = q*exp(1i*linspace(0,2*pi)');
confFill = [0.9,0.9,0.9];

zg = linspace(q,1,200)'.*exp(1i*linspace(0,2*pi,200));
nga = 10; ngb = 10;
zg(abs(zg-1)<1e-1)=nan;

for j = 5;%[1:3,5]
    if j ==1
    % Straining flow with circular wing
    [f,fp,a,d] = circularWing(q);
    vals = 4*[-linspace(0,1,20),linspace(0,1,20)];
    [potential,compVel] = strain(q,a);
    ax = [-4,4,-.2,4];
    elseif j ==2
    % Flat plate wing with circulatory flow
    aoa = pi/8;
    [f,fp,a,zt,d] = flatWing(aoa,q);
    vals = linspace(0,log(q),20)/2/pi;
    [potential,compVel] = circulatory;
    ax = [-.5,1.5,-.1,1];
    elseif j ==3
    % Zero aoa flat plate wing with motion
    [f,fp,a,zt,d] = flatWing(0,q);
    vals = linspace(-1,1,50);
    ni = 200;
    zInt = q*exp(1i*linspace(0,2*pi,ni)); zInt4 = permute(zInt,[1,4,3,2]);
    D = f(q); DDot = -1i; aDot = 0;
    I = -1/(2i*pi).*trapz(zInt,Mprime(zInt4,f,D,DDot,aDot)./zInt4,4);
    M = @(zVar) Mprime(zVar,f,D,DDot,aDot)+ I;
    movePot = @(zVar) zVar/pi.*trapz(zInt,M(zInt4).*PfunD(zVar./zInt4,q,N)./Pfun(zVar./zInt4,q,N)./zInt4.^2,4);
    potential  = movePot;
    ax = [-2,2,-.1,2];
    elseif j ==4
    % Circular arc wing with vortex
    gamma = .45*exp(.6i*pi/4);
    [f,fp,a,zt,d] = circularArcWing(gamma,q);
    [potential,compVel] = vortices(beta,q);
    vals = -linspace(0.01,.5,20);
    ax = .5*[-2,2,-.1,2];
    elseif j ==5
    % Centered circular arc with with uniform flow and Kutta condition   
    [f,fp,a,zt,d] = centeredCircularArcWing(phi,q);
    [uniPot,uniVel] = uniform(q,a);
    [cirPot,cirVel] = circulatory;
    circ = zt(2)*2i*pi*uniVel(zt(2));
    potential = @(zVar) circ*cirPot(zVar) + uniPot(zVar);
    vals = imag(potential(q))+linspace(-1,1,40);
    ax = [-1.5,2,-.1,1];
    end

% physical plot
f2 = figure(2);
clf
fill([-5 5 5 -5],[-.2 -.2 0 0],confFill,'EdgeColor',confFill);
hold on

contour(real(f(zg)),imag(f(zg)),imag(potential(zg)),vals,'k','LineWidth',.5)
if j ==5
    col = lines;
    v = [imag(potential(q)) imag(potential(q))];
    contour(real(f(zg)),imag(f(zg)),imag(potential(zg)),v,'LineColor','g','LineWidth',3)
end
if j == 4
   scatter(real(f(beta)),imag(f(beta)),150,'k','MarkerFaceColor','g') 
end
plot(f(c0),'b','LineWidth',7)
plot(f(c1),'r','LineWidth',7)
if j==1
fill(real(f(c1)),imag(f(c1)),confFill);
end

axis(ax);
hold off
axis off
ax = gca;
set(ax,'Units', 'normalized', 'Position', [0 0 1 1])
rat = range(ax.YLim)./range(ax.XLim);
f2.Position(3) = 1/rat*f2.Position(4);
%print(['z-strm',num2str(j)],'-dpng','-r100');

end

function Mp = Mprime(zVar,f,D,DDot,aDot)

Mp = imag(conj(DDot)*(f(zVar)-D) + 1i*aDot*abs(f(zVar)-D).^2);

end