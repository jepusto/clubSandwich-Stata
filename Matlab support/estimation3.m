clear
clc

load base_state_till_3
% legal = legal(1:27,:);
% beertaxa = beertaxa(1:27,:);
% mrate = mrate(1:27,:);
% pop = pop(1:27,:);
% state = state(1:27,:);

X = [legal beertaxa];

n = size(X,1);
%X = [ones(n,1) X];
p = size(X,2);

%pop = pop./100;
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
    %B_i = D_i'*(theta_i-X_i*M*X_i')*D_i;
    %B_i = D_i'*I_H_i*theta_i*D_i;
    
    %A_i = D_i'*sqrtm(pinv(B_i))*D_i;
    A_i = D_i'*sym_sqrt_MPinv(B_i)*D_i;
    XWAeeAWX = XWAeeAWX + X_i'*W_i*A_i*e_i*e_i'*A_i'*W_i*X_i;
    
    P(:,:,i) = I_H_i'*A_i*W_i*X_i*M;
    
    % X_i'*W_i*A_i
end

VCR = M*XWAeeAWX*M;
se = sqrt(diag(VCR))

%%
for k = 1:p
   c = zeros(p,1);
   c(k) = 1;
   
   top = 0;
   for  i = 1:size(clusters)
       top = top + (P(:,:,i)*c)'*theta*(P(:,:,i)*c);
   end
   
   bottom = 0;
   
   for  i = 1:size(clusters)
        for  j = 1:size(clusters)
            bottom = bottom + ((P(:,:,i)*c)'*theta*(P(:,:,j)*c))^2;
        end
   end
   
   df(k) = top^2/bottom;
    
end

df