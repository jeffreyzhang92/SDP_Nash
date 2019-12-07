function [e,x,y] = nash(A,B,num_its,alg)
%NASH Compute an approximate Nash equilibrium
% [e,x,y] = nash(A,B,[n],[alg])
% A is the payoff matrix of the row player.
% B is the payoff matrix of the column player.
% n is the number of iterations of the linearization algorithm.
% [alg] is the algorithm. Default is diagonal gap algorithm. 'SR' uses the
% square root algorithm.
% The epsilon returned is for a scaled version of the game which lies in
% [0,1].

if ~ismatrix(A) | ~ismatrix(B)
    error('A and B must be matrices.');
end
if ~isequal(size(A),size(B))
    error('A and B must be the same size.');
end
if ~exist('num_its','var')
    num_its=10;
end
if mod(num_its,1)~=0 | num_its < 0
    error('Number of iterations must be a nonnegative integer.');
end

if ~exist('alg','var') | ~strcmp(alg, 'SR')
    usesr = 0;
else
    usesr = 1;
end


if(max(max(A)) - min(min(A))~=0)
    A=(A - min(min(A)))/(max(max(A)) - min(min(A)));
end
if(max(max(B)) - min(min(B))~=0)
    B=(B - min(min(B)))/(max(max(B)) - min(min(B)));
end

n1 = size(A,1);
n2 = size(A,2);

cvx_begin sdp quiet

variable M(n1+n2,n1+n2) symmetric
variables x(n1) y(n2)

minimize(trace(M))

subject to

X = M(1:n1,1:n1);
Y = M(n1+1:n1+n2, n1+1:n1+n2);
Z = M(n1+1:n1+n2,1:n1);
P = M(1:n1,n1+1:n1+n2);
MXY = [M [x; y]; [x; y].' 1];
sum(x) == 1
sum(y) == 1
x >= 0
y >= 0
M(:) >= 0

MXY >= 0

sum(sum(X)) == 1
sum(sum(Y)) == 1
sum(sum(Z)) == 1

for(i = 1:n1)
    trace(A*Z) >= A(i,:)*y
    sum(X(:,i)) == sum(Z(:,i))
    sum(X(:,i)) == x(i)
end
for(i = 1:n2)
    trace(B*Z) >= x.' * B(:,i)
    sum(Y(i,:)) == sum(Z(i,:))
    sum(Y(i,:)) == y(i)
end

for i = 1:n1
    for j = 1:n1
        sum(A(i,:).*P(i,:)) >= sum(A(j,:).*P(i,:))-10e-7;
    end
end
for i = 1:n2
    for j = 1:n2
        sum(B(:,i).*P(:,i)) >= sum(B(:,j).*P(:,i))-10e-7;
    end
end

cvx_end


for k = 1:num_its
    
if usesr
    C = [diag(1./(diag(M+.0000001*eye(n1+n2)).^.5)) zeros(n1+n1,1);zeros(1,n1+n2),0];
else
    C = [eye(n1+n2) -[x;y];-[x;y].' 1];
end
    
cvx_begin sdp quiet

variable M2(n1+n2,n1+n2) symmetric
variables x2(n1) y2(n2)

minimize (trace(C*[M2 [x2; y2]; [x2; y2].' 1]))

subject to

X2 = M2(1:n1,1:n1);
Y2 = M2(n1+1:n1+n2, n1+1:n1+n2);
Z2 = M2(n1+1:n1+n2,1:n1);
P2 = M2(1:n1,n1+1:n1+n2);
MXY2 = [M2 [x2; y2]; [x2; y2].' 1];
sum(x2) == 1
sum(y2) == 1
x2 >= 0
y2 >= 0
M2(:) >= 0

MXY2 >= 0

sum(sum(X2)) == 1
sum(sum(Y2)) == 1
sum(sum(Z2)) == 1

for(i = 1:n1)
    trace(A*Z2) >= A(i,:)*y2
    sum(X2(:,i)) == sum(Z2(:,i))
    sum(X2(:,i)) == x2(i)
end
for(i = 1:n2)
    trace(B*Z2) >= x2.' * B(:,i)
    sum(Y2(i,:)) == sum(Z2(i,:))
    sum(Y2(i,:)) == y2(i)
end

for i = 1:n1
    for j = 1:n1
        sum(A(i,:).*P2(i,:)) >= sum(A(j,:).*P2(i,:))-10e-7;
    end
end
for i = 1:n2
    for j = 1:n2
        sum(B(:,i).*P2(:,i)) >= sum(B(:,j).*P2(:,i))-10e-7;
    end
end
cvx_end

u = max(A*y2)-A*y2;
v = max(B.'*x2)-B.'*x2;
vals(k,1) = cvx_optval;
vals(k,2) = (u.'*x2)^2 + (v.'*y2)^2;
vals(k,3) = (trace(A*Z2) - trace(y2*x2.'*A))/(max(max(A)) - min(min(A)));
vals(k,4) = (trace(B*Z2) - trace(y2*x2.'*B))/(max(max(B)) - min(min(B)));
M = M2;
x = x2;
y = y2;

end

e = max([A*y - x.'*A*y;B.'*x - x.'*B*y]);