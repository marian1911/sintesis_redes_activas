.ALIASES
X_U1A           U1A(+=N00407 -=N00733 V+=N00201 V-=N00215 OUT=N00733 ) CN
+@TP4_MONTECARLO.SCHEMATIC1(sch_1):INS32@OPAMP.LM324.Normal(chips)
X_U1B           U1B(+=N01000 -=OUT V+=N00201 V-=N00215 OUT=OUT ) CN
+@TP4_MONTECARLO.SCHEMATIC1(sch_1):INS73@OPAMP.LM324.Normal(chips)
V_V1            V1(+=N00201 -=0 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS123@SOURCE.VDC.Normal(chips)
V_V2            V2(+=0 -=N00215 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS139@SOURCE.VDC.Normal(chips)
R_R1            R1(1=N00800 2=N00407 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS282@ANALOG.R.Normal(chips)
R_R2            R2(1=N00415 2=N00800 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS298@ANALOG.R.Normal(chips)
R_R3            R3(1=0 2=N00415 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS314@ANALOG.R.Normal(chips)
R_R4            R4(1=N00415 2=IN ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS330@ANALOG.R.Normal(chips)
V_V3            V3(+=IN -=0 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS452@SOURCE.VSIN.Normal(chips)
C_C1            C1(1=0 2=N00407 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS564@ANALOG.C.Normal(chips)
C_C2            C2(1=N00800 2=N00733 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS771@ANALOG.C.Normal(chips)
C_C3            C3(1=N00733 2=N00996 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS927@ANALOG.C.Normal(chips)
C_C4            C4(1=N00996 2=N01000 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS943@ANALOG.C.Normal(chips)
R_R5            R5(1=0 2=N01000 ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS968@ANALOG.R.Normal(chips)
R_R6            R6(1=N00996 2=OUT ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS1184@ANALOG.R.Normal(chips)
R_R7            R7(1=0 2=OUT ) CN @TP4_MONTECARLO.SCHEMATIC1(sch_1):INS1241@ANALOG.R.Normal(chips)
_    _(in=IN)
_    _(out=OUT)
.ENDALIASES
