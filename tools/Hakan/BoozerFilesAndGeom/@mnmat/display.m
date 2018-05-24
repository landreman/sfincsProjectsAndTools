function display(A)
sz=size(A);
if length(sz)==2
  disp(['Forier size: 0<=m<=',num2str(size(A.c,1)-1),', ',...
        num2str(1-A.n0ind),'<=n<=',num2str(size(A.c,2)-A.n0ind)])
  disp(['Corresponding spatial discretisation size: ',num2str(sz(1)),'-by-',num2str(sz(2))])
  if strcmp(parity(A),'cos') || strcmp(parity(A),'both')
    disp('cos components = ')
    disp(A.c)
  end
  if strcmp(parity(A),'sin') || strcmp(parity(A),'both')
    disp('sin components = ')
    disp(A.s)
  end
else
  for ind=1:numel(A)
    tmp=A(:);
    tmp(ind)
  end
end 
 