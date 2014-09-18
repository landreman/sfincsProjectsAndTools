function varargout=fig(fignr,varargin)

sub=0;
strind=1;
firm=0;
for ind=1:nargin-1
  if ischar(varargin{ind})
    strs{strind}=varargin{ind};
    strind=strind+1;
  end
end
strind=strind-1;

if strind>=1
  firm1=strcmp(strs{1},'firm');
  sub1=strcmp(strs{1},'sub');
  nosub1=strcmp(strs{1},'nosub');
  if strind==2 
    firm2=strcmp(strs{1},'firm');
    sub2=strcmp(strs{1},'sub');
    nosub2=strcmp(strs{1},'nosub');
  else
    firm2=0;
    sub2=0; 
    nosub2=0; 
  end
  firm=firm1|firm2;
  sub=sub1|sub2;
  nosub=nosub1|nosub2;
  sub=sub & not(nosub);
end

renderer='painters';%'zbuffer';%'painters';%'OpenGL';%'zbuffer';


% dimensions becrux:0,43,1592,1078

scr_xwidth = 1592;  %1280 chalmers
scr_yheight=1078;  %1024 chalmers
scr_yskip=43;
yskip=325; % should by > scr_yskip
border=5;


if sub
  addthis=0;
  for ind=1:nargin-1
    if isnumeric(varargin{ind})
      addthis=varargin{ind};
    end
  end
  figure(500+addthis); %This number is hardcoded here
  fh=gcf;
  set(gcf,'Renderer',renderer)
  rows=3;
  cols=4;
  set(gcf,'Position',[0 0 scr_xwidth scr_yheight]);
  %col=mod(fignr-1,cols)+1;
  %row=(fignr-row)/rows;
  subplot(rows,cols,fignr)
  %set(gcf,'Position',[0 yskip scr_xwidth scr_yheight]);
else
  figure(fignr);
  fh=gcf;
  set(gcf,'Renderer',renderer)
  rows=2;
  cols=4;
  extrarow=1;
  
  pos=get(gcf,'Position');
  %pos=zeros(1,4);
  width=356; height=356;
  windowheight=height+84;
  windowwidth =width + 8;
  
  if (pos(3) ~= width) | firm
    %if firm
    
    pos(3)=width;
    pos(4)=height;
    
    %find which screen it is
    
    xscreen=0;%floor((pos(1)-border)/scr_xwidth);
    yscreen=0;%floor((pos(2)-border)/scr_yheight);
    xscreen=0;%scr_xwidth*xscreen;
    yscreen=0;%scr_yheight*yscreen;
    
    if fignr<=rows*cols
      pos(1)=xscreen+border+windowwidth*mod(fignr-1,cols);
      pos(2)=yscreen+yskip+windowheight*(mod(fignr-1,rows*cols)<cols);
    elseif and(extrarow,(fignr<=(rows+1)*cols))
      pos(1)=xscreen+border+windowwidth*mod(fignr-1,cols);
      pos(2)=yscreen+scr_yskip;    
      %  elseif fignr<=rows*cols*2 %At Chalmers it was possible to print to the
      %  adjacent screen
      %    pos(1)=xscreen+scr_xwidth+border+windowwidth*mod(fignr-1,cols);  
      %    pos(2)=yscreen+yskip+windowheight*(mod(fignr-1,rows*cols)<cols);
    else
      %pos(1)=xscreen+230;
      %pos(2)=yscreen+430;
    end
    figure(fignr);
    fh=gcf;
    set(gcf,'Position',[pos])
  else
    figure(fignr);
    fh=gcf;
  end  
end
if nargout>0
  varargout={fh};
end

set(gcf,'PaperPositionMode','auto','Color','w');%[0.8,0.8,0.8]);
colormap('jet')

