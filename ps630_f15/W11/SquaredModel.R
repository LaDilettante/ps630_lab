LDC$l1polity=LDC$l1polity+10

LDC$l1politysq=(LDC$l1polity)^2
LDC$l1fdisq=(LDC$l1fdi)^2
LDC$l1ushegsq=(LDC$l1usheg)^2

lm_newtar2=lm(newtar ~ l1polity + l1fdi + l1usheg + l1politysq + l1fdisq + l1ushegsq + factor(ctylabel)-1, data = LDC)
summary(lm_newtar2)

resettest(lm_newtar2, power = 2, data=LDC)
