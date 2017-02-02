clear
clc

load base_state_till_3

X = [legal beertaxa];

n = size(X,1);
p = size(X,2);

%% Choose 1:
%W = eye(n);
W = diag(pop);

%% Chosse 1:
%theta = eye(n);
theta = diag(1./pop);


%% beta and residuals:
Y = mrate;
b = inv(X'*W*X)*(X'*W*Y)
e = Y - X*b;

%% VCR:
%VCR = M*SUM(XWAeeAWX)*M

M = inv(X'*W*X);
I_H = eye(n)-X*M*X'*W;

% number of cluster:
clusters = unique(state);

XWAeeAWX = zeros(p,p);

for i = 1:size(clusters)
    c = clusters(i);
    
    loc = find(state==c);
    theta_i = theta(loc,loc);
    W_i = W(loc,loc);
    X_i = X(loc,:);
    e_i = e(loc,:);
    I_H_i = I_H(loc,:);
    D_i = chol(theta_i);
    
    B_i = D_i*I_H_i*theta*I_H_i'*D_i';
    
    A_i = D_i'*sqrtm(pinv(B_i))*D_i;
    XWAeeAWX = XWAeeAWX + X_i'*W_i*A_i*e_i*e_i'*A_i'*W_i*X_i;
end

VCR = M*XWAeeAWX*M;
se = sqrt(diag(VCR))