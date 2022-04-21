clear all
clc

%% Input

t = 0.1:0.1:4;
X = 0.1:0.01:1;

n = 2;
a = 0.1;
b = 0.2;

theta = 1;
mu = 0.7;

%% Jacobi

A = @(j) theta*(n-j).*((n-j-1)*a-1);
B = @(j) theta*(n-j).*((n-j-1)*b+mu);

P = [];
for i = 1:length(t)
    Q = [];
    for k = 0:n
        PB = prod(B(0:k-1));
        Sum = 0;
        for j = 0:k
            PA = prod(A(j)-A(0:j-1))*prod(A(j)-A(j+1:k));
            Sum = Sum + exp(t(i)*A(j))/PA;
        end
        Q = [Q Sum*PB];
    end
    P = [P; Q];
end

for i = 1:length(X)
    for j = 1:length(t)
        x(i,j) = X(i);
        y(i,j) = t(j);
        U(i,j) = sum(X(i).^(n:-1:0).*P(j,:));
    end
end

subplot(1,2,1)
surf(x,y,U)
subplot(1,2,2)
contourf(x,y,U,15)










