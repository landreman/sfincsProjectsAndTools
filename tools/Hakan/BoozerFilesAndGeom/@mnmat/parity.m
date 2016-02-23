function par=parity(Fmn)

hascos=any(abs(Fmn.c(find(not(isnan(Fmn.c)))))>10*eps);
hassin=any(abs(Fmn.s(find(not(isnan(Fmn.s)))))>10*eps);

if hascos && not(hassin)
  par='cos';
elseif not(hascos) && hassin
  par='sin';
else
  par='both';
end