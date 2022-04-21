clear all
clc

t = 0.1:0.1:1;
X = 0.1:0.1:1;

n = 2;
b = 0.5;
a = -b;

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


%% Monte Carlo

Npath = 10000;

dt = 0.001;
M = [];
for j = 1:length(t)
    T = 0:dt:t(j);
    R = repmat(X',1,Npath)';
    for i = 1:length(T)
        R = R + theta*(mu-R)*dt + sqrt(dt*2*b*theta*R.*(1-R)).*randn(Npath,length(X));
    end
    M = [M mean(R.^n)'];
	t(j)
end

C = abs(U-M);
[0 t; X' C]

subplot(1,2,1)
surf(x,y,C)
subplot(1,2,2)
contourf(x,y,C,15)



