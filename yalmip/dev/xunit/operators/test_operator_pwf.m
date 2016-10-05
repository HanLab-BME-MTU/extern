function test_operator_pwf

% Worked 20130322
sdpvar t
F1 = [-1 <= t <= 0];
F2 = [0 <= t <= 1];
f1 = (t+0.5)^2-1;
f2 = t^2;
z = pwf(f1,F1,f2,F2);
solvesdp([],z)
mbg_asserttolequal(double(t),-.5, 1e-4);

% Failed 20130322
sdpvar u
solvesdp([z <= u, -10 <= t <= 10],u)

% Failed 20130322
sdpvar s t
F1 = [-1 <= t <= 0];
G1 = [-1 <= s <= 0];
F2 = [0 <= t <= 1];
G2 = [0 <= s <= 1];
f1 = (t+0.5)^2-1;
g1 = (s+0.3)^2-1;
f2 = t^2;
g2 = s^2;
z = pwf(f1,F1,f2,F2);
w = pwf(g1,G1,g2,G2);
solvesdp([],z+w)

% Failed 20130322
sdpvar u v
solvesdp([z <= u, w <= v],u+v)

