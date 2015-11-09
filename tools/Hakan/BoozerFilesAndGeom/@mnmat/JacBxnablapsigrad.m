function out=JacBxnablapsigrad(arg,G,I)

%Assumes Boozer coordinates


gr=grad(arg);

out=gr(1)*G-gr(2)*I;