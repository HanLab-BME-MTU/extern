function test__misc_complexsdp1

i=sqrt(-1);
P = [4 1+2*i 3-i;1-2*i 3.5 0.8+2.3*i;3+i 0.8-2.3*i 4];
Zmanual = toeplitz([4.2827,0.8079+1.7342*sqrt(-1) 2.5574-0.7938*i])

Z = sdpvar(3,3,'toeplitz','complex');
Z = Z-sqrt(-1)*diag(imag(diag(Z)));

t = sdpvar(1,1);
e = Z(:)-P(:);
F = (Z >= 0);
F = F+ (norm(e) <= t);
sol = solvesdp(F,t)
assertTrue(sol.problem == 0);
assertTrue(norm(Zmanual-double(Z)) < 1e-4)

sol = solvesdp((Z >= 0),norm(e))
assertTrue(sol.problem == 0);
assertTrue(norm(Zmanual-double(Z)) < 1e-4)