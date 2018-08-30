function Fmn=removesin(Fmn)
Fmn.s(find(not(isnan(Fmn.s))))=0;