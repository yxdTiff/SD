\ENCODING=ISO-8859-1
\Problem name: baa99_pyomo

Minimize
 obj: 4 x1 + 2 x2 + 10 u1 + 10 u2 + 0.2 v1 + 0.2 v2 - 8 w11 - 4 w12 - 4 w22
Subject To
 c_e_d1_:              u1 + w11  = 106.6741630576
 c_e_d2_:              u2 + w12 + w22  = 102.63122844408
 c_e_s1_:              - x1 + v1 + w11 + w12  = 0
 c_e_s2_:              - x2 + v2 + w22  = 0
 c_e_ONE_VAR_CONSTANT: ONE_VAR_CONSTANT  = 1
Bounds
 0 <= x1 <= 217
 0 <= x2 <= 217
End
