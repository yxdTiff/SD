\ENCODING=ISO-8859-1
\Problem name: farmer

Minimize
 obj: 150 x1 + 230 x2 + 260 x3 + 238 y1 + 210 y2 - 170 w1 - 150 w2 - 36 w3
      - 10 w4
Subject To
 tLand:    x1 + x2 + x3 <= 500
 minWheat: 2.5 x1 + y1 - w1 >= 200
 minCorn:  3 x2 + y2 - w2 >= 240
 sSold:    - 20 x3 + w3 + w4 <= 0
 sFavor:   w3 <= 6000
End
