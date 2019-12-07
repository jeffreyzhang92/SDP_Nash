function [b] = nashbound(A,B,vec_x,vec_y,M_xy,M_xx,M_yy)
%NASHBOUND nashbound(A,B,vec_x,vec_y,M_xy,M_xx,M_yy) finds a lower bound for
%x^T(M_xy)y + x^T(M_xx)x + y^T(M_yy)y + x^T(vec_x) + y^T(vec_y)
%under any Nash equilibrium.
% If the argument matrix B is the string 'sym', then the game is assumed
% to be symmetric, and the bound is over symmetric equilibria

if ~ismatrix(A)
    error('A,B,M_xy,M_xx,M_yy must be matrices.');
end

n1 = size(A,1);
n2 = size(A,2);

if ischar(B) && strcmp(B,'sym')
    if n1 ~= n2
        error('A must be square in a symmetric game');
        end
    sym = 1;
    B = A';
else
    sym = 0;
end

if ~ismatrix(B)
    error('A,B,M_xy,M_xx,M_yy must be matrices.');
end
if ~isequal(size(A),size(B))
    error('A and B must be the same size.');
end

if ~exist('vec_x','var') | vec_x == 0
    vec_x = zeros(n1,1);
end
if ~exist('vec_y','var') | vec_y == 0
    vec_y = zeros(n2,1);
end
if ~exist('M_xx','var') | M_xx == 0
    M_xx = zeros(n1,n1);
end
if ~exist('M_xy','var') | M_xy == 0
    M_xy = zeros(n1,n2);
end
if ~exist('M_yy','var') | M_yy == 0
    M_yy = zeros(n2,n2);
end

if ~ismatrix(M_xy) | ~ismatrix(M_xx) | ~ismatrix(M_yy)
    error('A,B,M_xy,M_xx,M_yy must be matrices.');
end
if ~isvector(vec_x) | ~isvector(vec_y)
    error('vec_x and vec_y must be vectors.');
end

if(~isequal(size(M_xy),[n1 n2]))
    error('Dimension mismatch in M_xy.');
end
if(~isequal(size(M_xx),[n1 n1]))
    error('Dimension mismatch in M_xx.');
end
if(~isequal(size(M_yy),[n2 n2]))
    error('Dimension mismatch in M_yy.');
end
if(~isequal(size(vec_x),[n1 1]))
    error('Dimension mismatch in vec_x.');
end
if(~isequal(size(vec_y),[n2 1]))
    error('Dimension mismatch in vec_x.');
end

cvx_begin sdp quiet

variable M(n1+n2,n1+n2) symmetric
variables x(n1) y(n2)

minimize(trace(M*[M_xx M_xy/2; M_xy.'/2 M_yy])+x.'*vec_x+y.'*vec_y)
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
if sym == 1
    P == X;
    X == Y;
end

MXY >= 0

sum(sum(P)) == 1;
sum(P) == sum(Y);
sum(P,2) == sum(X,2);

trace(A*Z) >= A*sum(P)';
trace(B*Z) >= B'*sum(P,2);

cvx_end

b = cvx_optval;