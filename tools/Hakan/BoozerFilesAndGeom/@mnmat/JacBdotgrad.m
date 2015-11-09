function out=JacBdotgrad(arg,iota)

%Assumes Boozer coordinates

gr=grad(arg);

out=gr(1)*iota + gr(2);



