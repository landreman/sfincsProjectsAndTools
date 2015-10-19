function bool=iscommentline(str)

bool=0;

if str(1)=='c' || str(1)=='C' || str(1)=='?'
  bool=1;
elseif str(1)==' '
  if str(2)=='c' || str(2)=='C' || str(2)=='?'
    bool=1;
  end
end