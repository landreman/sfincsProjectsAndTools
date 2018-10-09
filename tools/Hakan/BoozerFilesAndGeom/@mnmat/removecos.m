function Fmn=removecos(Fmn)
Fmn.c(find(not(isnan(Fmn.c))))=0;