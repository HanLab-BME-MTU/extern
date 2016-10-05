function test_operator_optimizer9
sdpvar x y u z 

P1 = optimizer([x <= u,y <= z], -x-y,[],{u,z},[x;y]);
sol1 = P1{{2,3}};
mbg_asserttolequal(sol1,[2;3], 1e-4);

sol1 = P1{{[2 4],[3 7]}};
mbg_asserttolequal(sol1,[2 4;3 7], 1e-4);



