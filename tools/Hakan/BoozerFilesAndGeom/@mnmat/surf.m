function surf(in,varargin)%whichpart)

if nargin==2
  whichpart=varargin{1};
  if strcmp(whichpart,'c') || strcmp(whichpart,'cos')
    showcos=1;
  elseif strcmp(whichpart,'s') || strcmp(whichpart,'sin')
    showcos=0;
  else
    error('Unrecognised second input to surf.m for mnmat!')
  end
else
  nonnans_c=find(not(isnan(in.c)));
  nonnans_s=find(not(isnan(in.s)));
  if all(in.s(nonnans_s)==0)
    showcos=1;
  elseif all(in.c(nonnans_c)==0)
    showcos=0;
  else
    error(['Both cos and sin components exist!',...
           'You must choose which one to surf!']')
  end
end

mmax=in.m(end,1);
nmax=in.n(1,end);

%Matlab does not draw the last one so we have to add an extra dummy point.
c=[in.c,in.c(:,end)];
c=[c;c(end,:)];
s=[in.s,in.s(:,end)];
s=[s;s(end,:)];
m=[in.m,in.m(:,end)];
m=[m;m(end,:)+1];
n=[in.n,in.n(:,end)+1];
n=[n;n(end,:)];


if showcos
  c(find(abs(c)<10*eps))=0;
  surf(m-0.5,n-0.5,log10(abs(c)));
  axis([-0.5 mmax+0.5 -nmax-0.5 nmax+0.5])
  shading flat;
  view(0,90)
  colorbar
  title('lg |cosinus component|')
else
  s(find(abs(s)<10*eps))=0;
  surf(m-0.5,n-0.5,log10(abs(s)));
  axis([-0.5 mmax+0.5 -nmax-0.5 nmax+0.5])
  shading flat;
  view(0,90)
  colorbar
  title('lg |sinus component|')
end