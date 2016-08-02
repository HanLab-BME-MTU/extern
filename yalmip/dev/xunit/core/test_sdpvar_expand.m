function test_sdpvar_expand

yalmip('clear')
N = 1;
w = sdpvar(N, 1); %weight
previousPos = rand(N, 1);
%constraint = [0.0001 <= abs(w) <= 0.04, sum(abs(w)) == 1, abs(sum(w)) <= 0.04, sum(abs(w - previousPos)) <= 0.1];
constraint = [0.01 <= abs(w), abs(w) == 1, abs(w) == 1];
[~,~,~,m] = export(constraint,[],sdpsettings('solver',''));
assertTrue(length( m.binary_variables) == 1);

