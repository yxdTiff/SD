\ENCODING=ISO-8859-1
\Problem name: gscrc

Minimize
 obj: 1.5 F(P1,P2,Pr1,Truck,1)#0 + 1.5 F(P1,DC1,Pr1,Truck,1)#1
      + 1.5 F(P1,DC2,Pr1,Truck,1)#2 + 1.5 F(P2,P1,Pr1,Truck,1)#3
      + 1.5 F(P2,DC3,Pr1,Truck,1)#4 + 1.5 F(DC1,P1,Pr1,Truck,1)#5
      + 1.5 F(DC2,P1,Pr1,Truck,1)#6 + 1.5 F(DC2,DC3,Pr1,Truck,1)#7
      + 1.5 F(DC3,P2,Pr1,Truck,1)#8 + 1.5 F(DC3,DC2,Pr1,Truck,1)#9
      + 0.5 I(P1,Pr1,1)#10 + 0.5 I(P2,Pr1,1)#11 + 0.5 I(DC1,Pr1,1)#12
      + 0.5 I(DC2,Pr1,1)#13 + 0.5 I(DC3,Pr1,1)#14 + 2 S(P1,R1,Pr1,Truck,1)#15
      + 2 S(DC1,R1,Pr1,Truck,1)#16 + 2 S(DC1,R2,Pr1,Truck,1)#17
      + 2 S(DC1,R3,Pr1,Truck,1)#18 + 2 S(DC2,R4,Pr1,Truck,1)#19
      + 2 S(DC2,R5,Pr1,Truck,1)#20 + 2 S(DC3,R6,Pr1,Truck,1)#21
      + 500 SF(R1,Pr1,1)#24 + 500 SF(R2,Pr1,1)#25 + 500 SF(R3,Pr1,1)#26
      + 500 SF(R4,Pr1,1)#27 + 500 SF(R5,Pr1,1)#28 + 500 SF(R6,Pr1,1)#29
      + 1.5 F(P1,P2,Pr1,Truck,2)#30 + 1.5 F(P1,P2,Pr1,Truck,3)#31
      + F(P1,DC1,Pr1,Truck,2)#32 + F(P1,DC1,Pr1,Truck,3)#33
      + 1.5 F(P1,DC2,Pr1,Truck,2)#34 + 1.5 F(P1,DC2,Pr1,Truck,3)#35
      + 1.5 F(P2,P1,Pr1,Truck,2)#36 + 1.5 F(P2,P1,Pr1,Truck,3)#37
      + 1.5 F(P2,DC3,Pr1,Truck,2)#38 + 1.5 F(P2,DC3,Pr1,Truck,3)#39
      + 1.5 F(DC1,P1,Pr1,Truck,2)#40 + 1.5 F(DC1,P1,Pr1,Truck,3)#41
      + 1.5 F(DC2,P1,Pr1,Truck,2)#42 + 1.5 F(DC2,P1,Pr1,Truck,3)#43
      + 1.5 F(DC2,DC3,Pr1,Truck,2)#44 + 1.5 F(DC2,DC3,Pr1,Truck,3)#45
      + 1.5 F(DC3,P2,Pr1,Truck,2)#46 + 1.5 F(DC3,P2,Pr1,Truck,3)#47
      + 1.5 F(DC3,DC2,Pr1,Truck,2)#48 + 1.5 F(DC3,DC2,Pr1,Truck,3)#49
      + 0.5 I(P1,Pr1,2)#50 + 0.5 I(P1,Pr1,3)#52 + 0.5 I(P2,Pr1,2)#54
      + 0.5 I(P2,Pr1,3)#56 + 0.5 I(DC1,Pr1,2)#58 + 0.5 I(DC1,Pr1,3)#60
      + 0.5 I(DC2,Pr1,2)#62 + 0.5 I(DC2,Pr1,3)#64 + 0.5 I(DC3,Pr1,2)#66
      + 0.5 I(DC3,Pr1,3)#68 + 2 S(P1,R1,Pr1,Truck,2)#70
      + 2 S(P1,R1,Pr1,Truck,3)#71 + 2 S(DC1,R1,Pr1,Truck,2)#72
      + 2 S(DC1,R1,Pr1,Truck,3)#73 + 2 S(DC1,R2,Pr1,Truck,2)#74
      + 2 S(DC1,R2,Pr1,Truck,3)#75 + 2 S(DC1,R3,Pr1,Truck,2)#76
      + 2 S(DC1,R3,Pr1,Truck,3)#77 + 2 S(DC2,R4,Pr1,Truck,2)#78
      + 2 S(DC2,R4,Pr1,Truck,3)#79 + 2 S(DC2,R5,Pr1,Truck,2)#80
      + 2 S(DC2,R5,Pr1,Truck,3)#81 + 2 S(DC3,R6,Pr1,Truck,2)#82
      + 2 S(DC3,R6,Pr1,Truck,3)#83 + 500 SF(R1,Pr1,2)#92 + 500 SF(R1,Pr1,3)#93
      + 500 SF(R2,Pr1,2)#94 + 500 SF(R2,Pr1,3)#95 + 500 SF(R3,Pr1,2)#96
      + 500 SF(R3,Pr1,3)#97 + 500 SF(R4,Pr1,2)#98 + 500 SF(R4,Pr1,3)#99
      + 500 SF(R5,Pr1,2)#100 + 500 SF(R5,Pr1,3)#101 + 500 SF(R6,Pr1,2)#102
      + 500 SF(R6,Pr1,3)#103
Subject To
 massBalanceP1(Pr1,P1,1)#0:    F(P1,P2,Pr1,Truck,1)#0 + F(P1,DC1,Pr1,Truck,1)#1
                               + F(P1,DC2,Pr1,Truck,1)#2 + I(P1,Pr1,1)#10
                               + S(P1,R1,Pr1,Truck,1)#15 - W(P1,Pr1,1)#22  = 50
 massBalanceP1(Pr1,P2,1)#1:    F(P2,P1,Pr1,Truck,1)#3 + F(P2,DC3,Pr1,Truck,1)#4
                               + I(P2,Pr1,1)#11 - W(P2,Pr1,1)#23  = 50
 massBalanceDC1(Pr1,DC1,1)#2:  F(DC1,P1,Pr1,Truck,1)#5 + I(DC1,Pr1,1)#12
                               + S(DC1,R1,Pr1,Truck,1)#16
                               + S(DC1,R2,Pr1,Truck,1)#17
                               + S(DC1,R3,Pr1,Truck,1)#18  = 150
 massBalanceDC1(Pr1,DC2,1)#3:  F(DC2,P1,Pr1,Truck,1)#6
                               + F(DC2,DC3,Pr1,Truck,1)#7 + I(DC2,Pr1,1)#13
                               + S(DC2,R4,Pr1,Truck,1)#19
                               + S(DC2,R5,Pr1,Truck,1)#20  = 100
 massBalanceDC1(Pr1,DC3,1)#4:  F(DC3,P2,Pr1,Truck,1)#8
                               + F(DC3,DC2,Pr1,Truck,1)#9 + I(DC3,Pr1,1)#14
                               + S(DC3,R6,Pr1,Truck,1)#21  = 50
 capacity1(P1,Pr1,1)#5:        W(P1,Pr1,1)#22 <= 1000
 capacity1(P2,Pr1,1)#6:        W(P2,Pr1,1)#23 <= 800
 minimumInv1(P1,Pr1,1)#7:      I(P1,Pr1,1)#10 >= 20
 minimumInv1(P2,Pr1,1)#8:      I(P2,Pr1,1)#11 >= 10
 minimumInv1(DC1,Pr1,1)#9:     I(DC1,Pr1,1)#12 >= 150
 minimumInv1(DC2,Pr1,1)#10:    I(DC2,Pr1,1)#13 >= 100
 minimumInv1(DC3,Pr1,1)#11:    I(DC3,Pr1,1)#14 >= 50
 massBalanceC2(R1,Pr1,1)#12:   SF(R1,Pr1,1)#24  = 100
 massBalanceC2(R2,Pr1,1)#13:   SF(R2,Pr1,1)#25  = 100
 massBalanceC2(R3,Pr1,1)#14:   SF(R3,Pr1,1)#26  = 100
 massBalanceC2(R4,Pr1,1)#15:   SF(R4,Pr1,1)#27  = 100
 massBalanceC2(R5,Pr1,1)#16:   SF(R5,Pr1,1)#28  = 100
 massBalanceC2(R6,Pr1,1)#17:   SF(R6,Pr1,1)#29  = 100
 massBalanceP2(Pr1,P1,2)#18:   - F(P2,P1,Pr1,Truck,1)#3
                               - F(DC1,P1,Pr1,Truck,1)#5
                               - F(DC2,P1,Pr1,Truck,1)#6 - I(P1,Pr1,1)#10
                               + F(P1,P2,Pr1,Truck,2)#30
                               + F(P1,DC1,Pr1,Truck,2)#32
                               + F(P1,DC2,Pr1,Truck,2)#34 + I(P1,Pr1,2)#50
                               + S(P1,R1,Pr1,Truck,2)#70 - W(P1,Pr1,2)#84  = 0
 massBalanceP2(Pr1,P1,3)#19:   F(P1,P2,Pr1,Truck,3)#31
                               + F(P1,DC1,Pr1,Truck,3)#33
                               + F(P1,DC2,Pr1,Truck,3)#35
                               - F(P2,P1,Pr1,Truck,2)#36
                               - F(DC1,P1,Pr1,Truck,2)#40
                               - F(DC2,P1,Pr1,Truck,2)#42 - I(P1,Pr1,2)#50
                               + I(P1,Pr1,3)#52 + S(P1,R1,Pr1,Truck,3)#71
                               - W(P1,Pr1,3)#85  = 0
 massBalanceP2(Pr1,P2,2)#20:   - F(P1,P2,Pr1,Truck,1)#0
                               - F(DC3,P2,Pr1,Truck,1)#8 - I(P2,Pr1,1)#11
                               + F(P2,P1,Pr1,Truck,2)#36
                               + F(P2,DC3,Pr1,Truck,2)#38 + I(P2,Pr1,2)#54
                               - W(P2,Pr1,2)#86  = 0
 massBalanceP2(Pr1,P2,3)#21:   - F(P1,P2,Pr1,Truck,2)#30
                               + F(P2,P1,Pr1,Truck,3)#37
                               + F(P2,DC3,Pr1,Truck,3)#39
                               - F(DC3,P2,Pr1,Truck,2)#46 - I(P2,Pr1,2)#54
                               + I(P2,Pr1,3)#56 - W(P2,Pr1,3)#87  = 0
 massBalanceDC2(Pr1,DC1,2)#22: - F(P1,DC1,Pr1,Truck,1)#1 - I(DC1,Pr1,1)#12
                               + F(DC1,P1,Pr1,Truck,2)#40 + I(DC1,Pr1,2)#58
                               + S(DC1,R1,Pr1,Truck,2)#72
                               + S(DC1,R2,Pr1,Truck,2)#74
                               + S(DC1,R3,Pr1,Truck,2)#76  = 0
 massBalanceDC2(Pr1,DC1,3)#23: - F(P1,DC1,Pr1,Truck,2)#32
                               + F(DC1,P1,Pr1,Truck,3)#41 - I(DC1,Pr1,2)#58
                               + I(DC1,Pr1,3)#60 + S(DC1,R1,Pr1,Truck,3)#73
                               + S(DC1,R2,Pr1,Truck,3)#75
                               + S(DC1,R3,Pr1,Truck,3)#77  = 0
 massBalanceDC2(Pr1,DC2,2)#24: - F(P1,DC2,Pr1,Truck,1)#2
                               - F(DC3,DC2,Pr1,Truck,1)#9 - I(DC2,Pr1,1)#13
                               + F(DC2,P1,Pr1,Truck,2)#42
                               + F(DC2,DC3,Pr1,Truck,2)#44 + I(DC2,Pr1,2)#62
                               + S(DC2,R4,Pr1,Truck,2)#78
                               + S(DC2,R5,Pr1,Truck,2)#80  = 0
 massBalanceDC2(Pr1,DC2,3)#25: - F(P1,DC2,Pr1,Truck,2)#34
                               + F(DC2,P1,Pr1,Truck,3)#43
                               + F(DC2,DC3,Pr1,Truck,3)#45
                               - F(DC3,DC2,Pr1,Truck,2)#48 - I(DC2,Pr1,2)#62
                               + I(DC2,Pr1,3)#64 + S(DC2,R4,Pr1,Truck,3)#79
                               + S(DC2,R5,Pr1,Truck,3)#81  = 0
 massBalanceDC2(Pr1,DC3,2)#26: - F(P2,DC3,Pr1,Truck,1)#4
                               - F(DC2,DC3,Pr1,Truck,1)#7 - I(DC3,Pr1,1)#14
                               + F(DC3,P2,Pr1,Truck,2)#46
                               + F(DC3,DC2,Pr1,Truck,2)#48 + I(DC3,Pr1,2)#66
                               + S(DC3,R6,Pr1,Truck,2)#82  = 0
 massBalanceDC2(Pr1,DC3,3)#27: - F(P2,DC3,Pr1,Truck,2)#38
                               - F(DC2,DC3,Pr1,Truck,2)#44
                               + F(DC3,P2,Pr1,Truck,3)#47
                               + F(DC3,DC2,Pr1,Truck,3)#49 - I(DC3,Pr1,2)#66
                               + I(DC3,Pr1,3)#68 + S(DC3,R6,Pr1,Truck,3)#83
                                = 0
 massBalanceC2(R1,Pr1,2)#28:   S(P1,R1,Pr1,Truck,1)#15
                               + S(DC1,R1,Pr1,Truck,1)#16 + SF(R1,Pr1,2)#92
                               - SG(R1,Pr1,2)#104  = 100
 massBalanceC2(R1,Pr1,3)#29:   S(P1,R1,Pr1,Truck,2)#70
                               + S(DC1,R1,Pr1,Truck,2)#72 + SF(R1,Pr1,3)#93
                               - SG(R1,Pr1,3)#105  = 100
 massBalanceC2(R2,Pr1,2)#30:   S(DC1,R2,Pr1,Truck,1)#17 + SF(R2,Pr1,2)#94
                               - SG(R2,Pr1,2)#106  = 100
 massBalanceC2(R2,Pr1,3)#31:   S(DC1,R2,Pr1,Truck,2)#74 + SF(R2,Pr1,3)#95
                               - SG(R2,Pr1,3)#107  = 100
 massBalanceC2(R3,Pr1,2)#32:   S(DC1,R3,Pr1,Truck,1)#18 + SF(R3,Pr1,2)#96
                               - SG(R3,Pr1,2)#108  = 100
 massBalanceC2(R3,Pr1,3)#33:   S(DC1,R3,Pr1,Truck,2)#76 + SF(R3,Pr1,3)#97
                               - SG(R3,Pr1,3)#109  = 100
 massBalanceC2(R4,Pr1,2)#34:   S(DC2,R4,Pr1,Truck,1)#19 + SF(R4,Pr1,2)#98
                               - SG(R4,Pr1,2)#110  = 100
 massBalanceC2(R4,Pr1,3)#35:   S(DC2,R4,Pr1,Truck,2)#78 + SF(R4,Pr1,3)#99
                               - SG(R4,Pr1,3)#111  = 100
 massBalanceC2(R5,Pr1,2)#36:   S(DC2,R5,Pr1,Truck,1)#20 + SF(R5,Pr1,2)#100
                               - SG(R5,Pr1,2)#112  = 100
 massBalanceC2(R5,Pr1,3)#37:   S(DC2,R5,Pr1,Truck,2)#80 + SF(R5,Pr1,3)#101
                               - SG(R5,Pr1,3)#113  = 100
 massBalanceC2(R6,Pr1,2)#38:   S(DC3,R6,Pr1,Truck,1)#21 + SF(R6,Pr1,2)#102
                               - SG(R6,Pr1,2)#114  = 100
 massBalanceC2(R6,Pr1,3)#39:   S(DC3,R6,Pr1,Truck,2)#82 + SF(R6,Pr1,3)#103
                               - SG(R6,Pr1,3)#115  = 100
 capacity2(P1,Pr1,2)#40:       W(P1,Pr1,2)#84 + SW(P1,Pr1,2)#88  = 1000
 capacity2(P1,Pr1,3)#41:       W(P1,Pr1,3)#85 + SW(P1,Pr1,3)#89  = 1000
 capacity2(P2,Pr1,2)#42:       W(P2,Pr1,2)#86 + SW(P2,Pr1,2)#90  = 800
 capacity2(P2,Pr1,3)#43:       W(P2,Pr1,3)#87 + SW(P2,Pr1,3)#91  = 800
 minimumInv2(P1,Pr1,2)#44:     I(P1,Pr1,2)#50 - SI(P1,Pr1,2)#51  = 20
 minimumInv2(P1,Pr1,3)#45:     I(P1,Pr1,3)#52 - SI(P1,Pr1,3)#53  = 20
 minimumInv2(P2,Pr1,2)#46:     I(P2,Pr1,2)#54 - SI(P2,Pr1,2)#55  = 10
 minimumInv2(P2,Pr1,3)#47:     I(P2,Pr1,3)#56 - SI(P2,Pr1,3)#57  = 10
 minimumInv2(DC1,Pr1,2)#48:    I(DC1,Pr1,2)#58 - SI(DC1,Pr1,2)#59  = 150
 minimumInv2(DC1,Pr1,3)#49:    I(DC1,Pr1,3)#60 - SI(DC1,Pr1,3)#61  = 150
 minimumInv2(DC2,Pr1,2)#50:    I(DC2,Pr1,2)#62 - SI(DC2,Pr1,2)#63  = 100
 minimumInv2(DC2,Pr1,3)#51:    I(DC2,Pr1,3)#64 - SI(DC2,Pr1,3)#65  = 100
 minimumInv2(DC3,Pr1,2)#52:    I(DC3,Pr1,2)#66 - SI(DC3,Pr1,2)#67  = 50
 minimumInv2(DC3,Pr1,3)#53:    I(DC3,Pr1,3)#68 - SI(DC3,Pr1,3)#69  = 50
End
