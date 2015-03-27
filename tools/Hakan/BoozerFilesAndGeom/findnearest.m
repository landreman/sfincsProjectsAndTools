function inds=findnearest(vect,elems)


for y=1:length(elems)
  elem=elems(y);
  
  ind=find(vect==elem);
  if length(ind) ~= 1
    [m,ind]=min(abs(vect-elem));
    ind=ind(1);
  end
  inds(y)=ind;
end
