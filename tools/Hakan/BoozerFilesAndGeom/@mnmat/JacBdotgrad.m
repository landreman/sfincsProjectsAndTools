function out=JacBdotgrad(arg,iota,Nperiods)

%Assumes Boozer coordinates

gr=grad(arg,Nperiods);

out=gr(1)*iota + gr(2);



