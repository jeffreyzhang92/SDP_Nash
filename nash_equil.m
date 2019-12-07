m = 20;
n = 20;
A = random('unif',0,1, m, n);
B = random('unif',0,1, m, n);
M = [ zeros(m,m) A; B' zeros(n,n) ];

xi = ones(m+n,1);
eta = ones(m+n,1);
q = -ones(m+n,1);

maxiter = 10000;
sigma = q + M*xi - eta;
iter = 1;
while norm(sigma)>1e-8 && xi'*eta>1e-8 && iter<maxiter
    mu = 0.05*(xi'*eta)/(m+n);
    dxi = -(M + diag(eta./xi))\(M*xi + q - mu./xi); %solve
    deta = mu./xi - eta - (eta./xi).*dxi; %solve
    theta = 1/max(-[dxi; deta]./[xi; eta]); %max possible step size
    theta = 0.95*theta; %plus some wiggle room
    xi = xi + theta*dxi;
    eta = eta + theta*deta;
    sigma = q + M*xi - eta;
    iter = iter+1;
end
y = xi(1:m);
x = xi(m+1:end);
ystar = y/sum(y);
xstar = x/sum(x);

if iter == maxiter 
   sprintf('iteration limit. infeas=%s, noncomplimentarity=%s',norm(sigma), xi'*eta)
else
    A
    B
    xstar
    ystar
    sprintf('iterations = %d', iter)
    sprintf('y*T A x* = %f <= %f = min(A x*)', ystar'*A*xstar, min(A*xstar))
    sprintf('y*T A x* = %f <= %f = min(y*T B)', ystar'*B*xstar, min(ystar'*B))
end

