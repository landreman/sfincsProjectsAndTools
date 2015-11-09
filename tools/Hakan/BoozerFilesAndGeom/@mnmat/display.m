function display(A)
sz=size(A);
if length(sz)==2
  disp(['size ',num2str(sz(1)),'-by-',num2str(sz(2))])
  disp('cos components = ')
  disp(A.c)
  disp('sin components = ')
  disp(A.s)
else
  for ind=1:numel(A)
    tmp=A(:);
    tmp(ind)
  end
end 
 