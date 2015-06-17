\ENCODING=ISO-8859-1
\Problem name: gsc.cor

Minimize
 totalCost: 1.5 F(P1,P2,Pr1,Truck,1)#0 + 1.5 F(P1,DC1,Pr1,Truck,1)#1
            + 1.5 F(P1,DC2,Pr1,Truck,1)#2 + 1.5 F(P2,P1,Pr1,Truck,1)#3
            + 1.5 F(P2,DC3,Pr1,Truck,1)#4 + 1.5 F(DC1,P1,Pr1,Truck,1)#5
            + 1.5 F(DC2,P1,Pr1,Truck,1)#6 + 1.5 F(DC2,DC3,Pr1,Truck,1)#7
            + 1.5 F(DC3,P2,Pr1,Truck,1)#8 + 1.5 F(DC3,DC2,Pr1,Truck,1)#9
            + 0.5 I(P1,Pr1,1)#10 + 0.5 I(P2,Pr1,1)#11 + 0.5 I(DC1,Pr1,1)#12
            + 0.5 I(DC2,Pr1,1)#13 + 0.5 I(DC3,Pr1,1)#14
            + 2 S(P1,R1,Pr1,Truck,1)#15 + 2 S(DC1,R1,Pr1,Truck,1)#16
            + 2 S(DC1,R2,Pr1,Truck,1)#17 + 2 S(DC1,R3,Pr1,Truck,1)#18
            + 2 S(DC2,R4,Pr1,Truck,1)#19 + 2 S(DC2,R5,Pr1,Truck,1)#20
            + 2 S(DC3,R6,Pr1,Truck,1)#21 + 500 SF(R1,Pr1,1)#24
            + 500 SF(R2,Pr1,1)#25 + 500 SF(R3,Pr1,1)#26 + 500 SF(R4,Pr1,1)#27
            + 500 SF(R5,Pr1,1)#28 + 500 SF(R6,Pr1,1)#29
            + 1.5 F(P1,P2,Pr1,Truck,2)#30 + 1.5 F(P1,P2,Pr1,Truck,3)#31
            + 1.5 F(P1,DC1,Pr1,Truck,2)#32 + 1.5 F(P1,DC1,Pr1,Truck,3)#33
            + 1.5 F(P1,DC2,Pr1,Truck,2)#34 + 1.5 F(P1,DC2,Pr1,Truck,3)#35
            + 1.5 F(P2,P1,Pr1,Truck,2)#36 + 1.5 F(P2,P1,Pr1,Truck,3)#37
            + 1.5 F(P2,DC3,Pr1,Truck,2)#38 + 1.5 F(P2,DC3,Pr1,Truck,3)#39
            + 1.5 F(DC1,P1,Pr1,Truck,2)#40 + 1.5 F(DC1,P1,Pr1,Truck,3)#41
            + 1.5 F(DC2,P1,Pr1,Truck,2)#42 + 1.5 F(DC2,P1,Pr1,Truck,3)#43
            + 1.5 F(DC2,DC3,Pr1,Truck,2)#44 + 1.5 F(DC2,DC3,Pr1,Truck,3)#45
            + 1.5 F(DC3,P2,Pr1,Truck,2)#46 + 1.5 F(DC3,P2,Pr1,Truck,3)#47
            + 1.5 F(DC3,DC2,Pr1,Truck,2)#48 + 1.5 F(DC3,DC2,Pr1,Truck,3)#49
            + 0.5 I(P1,Pr1,2)#50 + 0.5 I(P1,Pr1,3)#51 + 0.5 I(P2,Pr1,2)#52
            + 0.5 I(P2,Pr1,3)#53 + 0.5 I(DC1,Pr1,2)#54 + 0.5 I(DC1,Pr1,3)#55
            + 0.5 I(DC2,Pr1,2)#56 + 0.5 I(DC2,Pr1,3)#57 + 0.5 I(DC3,Pr1,2)#58
            + 0.5 I(DC3,Pr1,3)#59 + 2 S(P1,R1,Pr1,Truck,2)#60
            + 2 S(P1,R1,Pr1,Truck,3)#61 + 2 S(DC1,R1,Pr1,Truck,2)#62
            + 2 S(DC1,R1,Pr1,Truck,3)#63 + 2 S(DC1,R2,Pr1,Truck,2)#64
            + 2 S(DC1,R2,Pr1,Truck,3)#65 + 2 S(DC1,R3,Pr1,Truck,2)#66
            + 2 S(DC1,R3,Pr1,Truck,3)#67 + 2 S(DC2,R4,Pr1,Truck,2)#68
            + 2 S(DC2,R4,Pr1,Truck,3)#69 + 2 S(DC2,R5,Pr1,Truck,2)#70
            + 2 S(DC2,R5,Pr1,Truck,3)#71 + 2 S(DC3,R6,Pr1,Truck,2)#72
            + 2 S(DC3,R6,Pr1,Truck,3)#73 + 500 SF(R1,Pr1,2)#78
            + 500 SF(R1,Pr1,3)#79 + 500 SF(R2,Pr1,2)#80 + 500 SF(R2,Pr1,3)#81
            + 500 SF(R3,Pr1,2)#82 + 500 SF(R3,Pr1,3)#83 + 500 SF(R4,Pr1,2)#84
            + 500 SF(R4,Pr1,3)#85 + 500 SF(R5,Pr1,2)#86 + 500 SF(R5,Pr1,3)#87
            + 500 SF(R6,Pr1,2)#88 + 500 SF(R6,Pr1,3)#89
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
 massBalanceP2(Pr1,P1,2)#12:   - F(P2,P1,Pr1,Truck,1)#3
                               - F(DC1,P1,Pr1,Truck,1)#5
                               - F(DC2,P1,Pr1,Truck,1)#6 - I(P1,Pr1,1)#10
                               + F(P1,P2,Pr1,Truck,2)#30
                               + F(P1,DC1,Pr1,Truck,2)#32
                               + F(P1,DC2,Pr1,Truck,2)#34 + I(P1,Pr1,2)#50
                               + S(P1,R1,Pr1,Truck,2)#60 - W(P1,Pr1,2)#74  = 0
 massBalanceP2(Pr1,P1,3)#13:   F(P1,P2,Pr1,Truck,3)#31
                               + F(P1,DC1,Pr1,Truck,3)#33
                               + F(P1,DC2,Pr1,Truck,3)#35
                               - F(P2,P1,Pr1,Truck,2)#36
                               - F(DC1,P1,Pr1,Truck,2)#40
                               - F(DC2,P1,Pr1,Truck,2)#42 - I(P1,Pr1,2)#50
                               + I(P1,Pr1,3)#51 + S(P1,R1,Pr1,Truck,3)#61
                               - W(P1,Pr1,3)#75  = 0
 massBalanceP2(Pr1,P2,2)#14:   - F(P1,P2,Pr1,Truck,1)#0
                               - F(DC3,P2,Pr1,Truck,1)#8 - I(P2,Pr1,1)#11
                               + F(P2,P1,Pr1,Truck,2)#36
                               + F(P2,DC3,Pr1,Truck,2)#38 + I(P2,Pr1,2)#52
                               - W(P2,Pr1,2)#76  = 0
 massBalanceP2(Pr1,P2,3)#15:   - F(P1,P2,Pr1,Truck,2)#30
                               + F(P2,P1,Pr1,Truck,3)#37
                               + F(P2,DC3,Pr1,Truck,3)#39
                               - F(DC3,P2,Pr1,Truck,2)#46 - I(P2,Pr1,2)#52
                               + I(P2,Pr1,3)#53 - W(P2,Pr1,3)#77  = 0
 massBalanceDC2(Pr1,DC1,2)#16: - F(P1,DC1,Pr1,Truck,1)#1 - I(DC1,Pr1,1)#12
                               + F(DC1,P1,Pr1,Truck,2)#40 + I(DC1,Pr1,2)#54
                               + S(DC1,R1,Pr1,Truck,2)#62
                               + S(DC1,R2,Pr1,Truck,2)#64
                               + S(DC1,R3,Pr1,Truck,2)#66  = 0
 massBalanceDC2(Pr1,DC1,3)#17: - F(P1,DC1,Pr1,Truck,2)#32
                               + F(DC1,P1,Pr1,Truck,3)#41 - I(DC1,Pr1,2)#54
                               + I(DC1,Pr1,3)#55 + S(DC1,R1,Pr1,Truck,3)#63
                               + S(DC1,R2,Pr1,Truck,3)#65
                               + S(DC1,R3,Pr1,Truck,3)#67  = 0
 massBalanceDC2(Pr1,DC2,2)#18: - F(P1,DC2,Pr1,Truck,1)#2
                               - F(DC3,DC2,Pr1,Truck,1)#9 - I(DC2,Pr1,1)#13
                               + F(DC2,P1,Pr1,Truck,2)#42
                               + F(DC2,DC3,Pr1,Truck,2)#44 + I(DC2,Pr1,2)#56
                               + S(DC2,R4,Pr1,Truck,2)#68
                               + S(DC2,R5,Pr1,Truck,2)#70  = 0
 massBalanceDC2(Pr1,DC2,3)#19: - F(P1,DC2,Pr1,Truck,2)#34
                               + F(DC2,P1,Pr1,Truck,3)#43
                               + F(DC2,DC3,Pr1,Truck,3)#45
                               - F(DC3,DC2,Pr1,Truck,2)#48 - I(DC2,Pr1,2)#56
                               + I(DC2,Pr1,3)#57 + S(DC2,R4,Pr1,Truck,3)#69
                               + S(DC2,R5,Pr1,Truck,3)#71  = 0
 massBalanceDC2(Pr1,DC3,2)#20: - F(P2,DC3,Pr1,Truck,1)#4
                               - F(DC2,DC3,Pr1,Truck,1)#7 - I(DC3,Pr1,1)#14
                               + F(DC3,P2,Pr1,Truck,2)#46
                               + F(DC3,DC2,Pr1,Truck,2)#48 + I(DC3,Pr1,2)#58
                               + S(DC3,R6,Pr1,Truck,2)#72  = 0
 massBalanceDC2(Pr1,DC3,3)#21: - F(P2,DC3,Pr1,Truck,2)#38
                               - F(DC2,DC3,Pr1,Truck,2)#44
                               + F(DC3,P2,Pr1,Truck,3)#47
                               + F(DC3,DC2,Pr1,Truck,3)#49 - I(DC3,Pr1,2)#58
                               + I(DC3,Pr1,3)#59 + S(DC3,R6,Pr1,Truck,3)#73
                                = 0
 massBalanceC2(R1,Pr1,1)#22:   SF(R1,Pr1,1)#24  = 100
 massBalanceC2(R1,Pr1,2)#23:   S(P1,R1,Pr1,Truck,1)#15
                               + S(DC1,R1,Pr1,Truck,1)#16 + SF(R1,Pr1,2)#78
                                = 110
 massBalanceC2(R1,Pr1,3)#24:   S(P1,R1,Pr1,Truck,2)#60
                               + S(DC1,R1,Pr1,Truck,2)#62 + SF(R1,Pr1,3)#79
                                = 100
 massBalanceC2(R2,Pr1,1)#25:   SF(R2,Pr1,1)#25  = 100
 massBalanceC2(R2,Pr1,2)#26:   S(DC1,R2,Pr1,Truck,1)#17 + SF(R2,Pr1,2)#80
                                = 110
 massBalanceC2(R2,Pr1,3)#27:   S(DC1,R2,Pr1,Truck,2)#64 + SF(R2,Pr1,3)#81
                                = 105
 massBalanceC2(R3,Pr1,1)#28:   SF(R3,Pr1,1)#26  = 100
 massBalanceC2(R3,Pr1,2)#29:   S(DC1,R3,Pr1,Truck,1)#18 + SF(R3,Pr1,2)#82
                                = 100
 massBalanceC2(R3,Pr1,3)#30:   S(DC1,R3,Pr1,Truck,2)#66 + SF(R3,Pr1,3)#83
                                = 110
 massBalanceC2(R4,Pr1,1)#31:   SF(R4,Pr1,1)#27  = 100
 massBalanceC2(R4,Pr1,2)#32:   S(DC2,R4,Pr1,Truck,1)#19 + SF(R4,Pr1,2)#84
                                = 110
 massBalanceC2(R4,Pr1,3)#33:   S(DC2,R4,Pr1,Truck,2)#68 + SF(R4,Pr1,3)#85
                                = 100
 massBalanceC2(R5,Pr1,1)#34:   SF(R5,Pr1,1)#28  = 100
 massBalanceC2(R5,Pr1,2)#35:   S(DC2,R5,Pr1,Truck,1)#20 + SF(R5,Pr1,2)#86
                                = 100
 massBalanceC2(R5,Pr1,3)#36:   S(DC2,R5,Pr1,Truck,2)#70 + SF(R5,Pr1,3)#87
                                = 110
 massBalanceC2(R6,Pr1,1)#37:   SF(R6,Pr1,1)#29  = 100
 massBalanceC2(R6,Pr1,2)#38:   S(DC3,R6,Pr1,Truck,1)#21 + SF(R6,Pr1,2)#88
                                = 110
 massBalanceC2(R6,Pr1,3)#39:   S(DC3,R6,Pr1,Truck,2)#72 + SF(R6,Pr1,3)#89
                                = 100
 capacity2(P1,Pr1,2)#40:       W(P1,Pr1,2)#74 <= 1000
 capacity2(P1,Pr1,3)#41:       W(P1,Pr1,3)#75 <= 1000
 capacity2(P2,Pr1,2)#42:       W(P2,Pr1,2)#76 <= 800
 capacity2(P2,Pr1,3)#43:       W(P2,Pr1,3)#77 <= 800
 minimumInv2(P1,Pr1,2)#44:     I(P1,Pr1,2)#50 >= 20
 minimumInv2(P1,Pr1,3)#45:     I(P1,Pr1,3)#51 >= 20
 minimumInv2(P2,Pr1,2)#46:     I(P2,Pr1,2)#52 >= 10
 minimumInv2(P2,Pr1,3)#47:     I(P2,Pr1,3)#53 >= 10
 minimumInv2(DC1,Pr1,2)#48:    I(DC1,Pr1,2)#54 >= 150
 minimumInv2(DC1,Pr1,3)#49:    I(DC1,Pr1,3)#55 >= 150
 minimumInv2(DC2,Pr1,2)#50:    I(DC2,Pr1,2)#56 >= 100
 minimumInv2(DC2,Pr1,3)#51:    I(DC2,Pr1,3)#57 >= 100
 minimumInv2(DC3,Pr1,2)#52:    I(DC3,Pr1,2)#58 >= 50
 minimumInv2(DC3,Pr1,3)#53:    I(DC3,Pr1,3)#59 >= 50
End
