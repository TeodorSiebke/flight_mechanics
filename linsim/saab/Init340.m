%Init340.m 340
% ADDB ref.
% xref = 25%C = STA 14.884 m
% zref = 0.254 m below centreline = WL 2.286 m 
%
% xcg = (xcgmac-25)*2.407e-2 + 14.884
% xcgmac=25
q = QA(1) ;
m = MASS(1) ;
Iy = IY(1) ;
u0=TAS(1);
theta0=THEG(1)*pi/180;

ALONG = [ GV11(1) GV12(1) GV13(1) GV14(1)
          GV21(1) GV22(1) GV23(1) GV24(1)
          GV31(1) GV32(1) GV33(1) GV34(1)
          GV41(1) GV42(1) GV43(1) GV44(1)] ;

BLONG = [ BV11(1) BV12(1)
          BV21(1) BV22(1)
          BV31(1) BV32(1)
          BV41(1) BV42(1)  ] ;

ALATE = [ GH11(1) GH12(1) GH13(1) GH14(1)
          GH21(1) GH22(1) GH23(1) GH24(1)
          GH31(1) GH32(1) GH33(1) GH34(1)
          GH41(1) GH42(1) GH43(1) GH44(1)] ;

BLATE = [ BH11(1) BH12(1)
          BH21(1) BH22(1)
          BH31(1) BH32(1)
          BH41(1) BH42(1)  ] ;
