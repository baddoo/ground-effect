% This function checks that the residues are correct

addpath('../maps','../primeP','../complex-potentials')

q = .8*rand; phi = pi*rand;
alpha = 2*pi*rand;
gamma = (q+(1-q)*rand)*exp(2i*pi*rand);

for j = 1:5
    switch j
        case 1
            [f,fd,a,d] = circularWing(q);
        case 2
            [f,fd,a,zt,d] = flatWing(alpha,q);
        case 3
            [f,fd,a,zt,d] = flatWing(0,q);
        case 4
            [f,fd,a,zt,d] = centeredCircularArcWing(phi,q);
        case 5
            [f,fd,a,zt,d] = circularArcWing(gamma,q);
    end

% Define line limiting to zeta = 1.
zLim = linspace(q,1); zLim(end) = [];
approx = a./(zLim-1);
err = f(zLim) - approx; 

% Plot error -- should go to a constant.
plot(abs(err),'LineWidth',3);
hold on

end

hold off
