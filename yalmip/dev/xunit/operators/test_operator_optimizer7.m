function test_operator_optimizer7

X = sdpvar(2,4,2);
Y = sdpvar(2,4,2);
P = optimizer([],sum(sum(abs(X(:)-Y(:)))),[],Y,X);
Z = reshape(magic(4),2,4,2);
U = P{Z};
mbg_asserttrue(norm(Z(:)-U(:))<1e-7);


% Test nD parameter in cells
X = sdpvar(2,4,2);
Y1 = sdpvar(2,4,2);
Y2 = sdpvar(2,4,2);
P = optimizer([],sum(sum(abs(X(:)-Y1(:))))+sum(sum(abs(X(:)-Y2(:)))),[],{Y1,Y2},X);
Z = reshape(magic(4),2,4,2);
U = P{{Z,Z}};
mbg_asserttrue( norm(Z(:)-U(:))<1e-8);

% Test nD outputs in cells
X = sdpvar(2,4,2);
Y1 = sdpvar(2,4,2);
P = optimizer([],sum(sum(abs(X(:)-Y1(:)))),[],Y1,{X,2*X});
Z = reshape(magic(4),2,4,2);
U = P{Z};
mbg_asserttrue( norm(Z(:)-U{1}(:))<1e-7);
mbg_asserttrue( norm(2*Z(:)-U{2}(:))<1e-7);

