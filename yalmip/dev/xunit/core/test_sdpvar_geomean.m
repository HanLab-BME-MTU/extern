function test_sdpvar_geomean

% Test real vector geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
A = randn(15,2);
b = rand(15,1)*5;
x = sdpvar(2,1);
obj = geomean(b-A*x);
solvesdp([],-obj);
assertElementsAlmostEqual(double(x'), [-0.05519469470525   0.26970610928222],'absolute', 1e-3);
assertElementsAlmostEqual(double(obj), 1.83896843735621,'absolute', 1e-3);

% Test real vector geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
A = randn(16,2);
b = rand(16,1)*5;
x = sdpvar(2,1);
obj = geomean(b-A*x);
solvesdp([],-obj);
assertElementsAlmostEqual(double(x'), [ -0.01148934254297  -0.20720944929269],'absolute', 1e-3);
assertElementsAlmostEqual(double(obj), 1.93924577959868,'absolute', 1e-3);

% Test real vector geomean, length == 1
randn('seed',1234);
rand('seed',1234);
A = randn(1,2);
b = rand(1,1)*5;
x = sdpvar(2,1);
obj = geomean(b-A*x);
sol = solvesdp([],-obj);
assertTrue(sol.problem==2);

% Test real matrix geomean, length ~2^n
randn('seed',1234);
rand('seed',1234);
D = randn(5,5);
P = sdpvar(5,5);
obj = geomean(P);
solvesdp((P <= D*D'),-obj);
assertElementsAlmostEqual(double(obj), 2.00333629658259, 'absolute',1e-4);
 
% Test real matrix geomean, length == 2^n
randn('seed',1234);
rand('seed',1234);
D = randn(8,8);
P = sdpvar(8,8);
obj = geomean(P);
solvesdp((P <= D*D'),-obj);
assertElementsAlmostEqual(double(obj), 3.32199302165511,'absolute', 1e-4);

% Test real matrix geomean, length == 2
randn('seed',1234);
rand('seed',1234);
D = randn(2,2);
P = sdpvar(2,2);
obj = geomean(P);
solvesdp((P <= D*D'),-obj);
assertElementsAlmostEqual(double(obj),  2.02896175488410,'absolute',1e-4);