function dk=read_dkes_dkfile(filename,varargin)


if isstruct(filename) %a dk struct can also be given direct, if only the
                      %interpolation is preparation is wanted, see below
  dk=filename;
  if nargin>1
    dk.Einterp_type=varargin{1};
  else
    dk.Einterp_type='linear';
  end

else
  if nargin>1
    dk.Einterp_type=varargin{1};
  else
    dk.Einterp_type='linear';
  end

  if filename(end-2:end)~='.dk'
    filename=[filename,'.dk'] %also accept only the nickname
  end
  %dk.filename=filename;
  
  fid = fopen(filename);
  if fid<0
    fid
    error(['Unable to open file:  ',filename])
  end

  %I assume that there is at least one comment line
  linestr=fgetl(fid);
  commentline=1;
  dk.comments{commentline}=linestr;
  while iscommentline(linestr)
    commentline=commentline+1;
    linestr=fgetl(fid);
    dk.comments{commentline}=linestr;
  end

  rind=0;
  eof=0;

  while not(eof)
    rind=rind+1;
    header1=sscanf(linestr,'%f %f %f %f %f %f %f',7);
    foundmaindata=0;
    while not(foundmaindata)
      linestr=fgetl(fid);
      foundmaindata=isempty(find(linestr=='c'));
      if not(foundmaindata)
        foundcfit=0;
        for ind=1:length(linestr)-3
          if linestr(ind:ind+3)=='cfit'
            foundcfit=1;
          end
        end
        if foundcfit
          header2=sscanf(linestr(6:end),'%f %f %f %f %f',7);
          dk.cfit.eps_eff(rind)  =header2(1);
          dk.cfit.g11_0(rind)    =header2(2);
          dk.cfit.efield_u(rind) =header2(3);
          dk.cfit.g11_sq(rind)   =header2(4);
          dk.cfit.ex_er(rind)    =header2(5);
        end
      end
    end
    %if rind==1
    %  fgetl(fid);  %Skip variable name line only shown in beginning of file
    %               % c         eps_eff     g11_ft      efield_u    g11_er    ex_er
    %end
    %linestr=fgetl(fid);
    %header2=sscanf(linestr(6:end),'%f %f %f %f %f',7);

    %---------------------------------------------------------------------------
    %HS,          HM notation
    %         "standard" input:
    %r            xr     minor radius in m
    %R            xrm    major radius in m (used for DKES computations)
    %B0           b0     magnetic field in T (used for DKES computations)
    %iota         iota   rotational transform
    %          additional (new) input:
    %NormTorCurv  xkn    norm. tor. curvature 
    %                    (for extrapolation of g11 in PS-regime)
    %FracTrap     ftp    trapped part. fraction 
    %B2av         b2av   flux-surface average of B^2
    %
    % HM also defines
    %             r      minor radius in cm
    %             rm     major radius in cm
    %
    %---------------------------------------------------------------------------

    dk.r(rind)          =header1(1);  %HM: arr(rind) 
    dk.R(rind)          =header1(2);  %HM: xrm
    dk.B0(rind)         =header1(3);  %HM: b0
    dk.iota(rind)       =header1(4);  %HM: air(rind)
    dk.NormTorCurv(rind)=header1(5);  %HM: akn(rind)
    dk.FracTrap(rind)   =header1(6);  %HM: aftp(rind)
    if length(header1)>6
      dk.B2av(rind)       =header1(7);  %HM: ab2av(rind)
    else
      dk.B2av(rind)=NaN;
    end
    
    dk.B02(rind)=dk.B0(rind)^2;

    % normalize efield to "resonance" value
    fnrm_e  = 1/(dk.B0(rind)*dk.r(rind));
    
    
    recind=0;
    
    if length(linestr)>=2
      while linestr(2)=='c' %skip the c line
        linestr=fgetl(fid);  %readnext line
      end        
    end

    while linestr(1)~='e' && linestr(2)~='e'
      if length(linestr)>=2
        while linestr(2)=='c' %skip the c line
          linestr=fgetl(fid);  %readnext line
        end        
      end
      if linestr(1)=='e' || linestr(2)=='e'
        error('unforeseen!')
        %Todo: Check what happens if the last entry
        %before 'e' is commented out!
      end
      recind=recind+1;
      line1=sscanf(linestr,'%f %f %f %f %f',7);
      linestr=fgetl(fid);
      line2=sscanf(linestr(5:end),'%d %d %f %f %f',7);
      
      % HM comment:
      %     read transport coefficient and normalize with respect to b0
      %    
      %     DKES logics:   (corrected in April 2008)
      %     g11 propto 1/b0**2, g13 independent of b0, g33 propto b0**2
      %     renormalisation:  in entry tc_dkes_0
      %     new definition of D33:   (February 2009)
      %     "old" D33 is replaced by D33/b2av  (see also dkes3.f)
      %     new:  1.5*cmul*d33 -> 1  for  large cmul

      dk.data{rind}.cmul(recind)   =line1(1); %nu/v %HM: xhnu
      dk.data{rind}.EovervB(recind)=line1(2); %E/vB %HM: xref
      dk.data{rind}.g11_i(recind)  =line1(3); %mean of bounds
      dk.data{rind}.g13_i(recind)  =line1(4); %mean of bounds
      dk.data{rind}.g33_i(recind)  =line1(5); %mean of bounds

      %dk.data{rind}.axnu(recind)   = ...
      %    log10(dk.data{rind}.xhnu(recind));
      %xref_l=1.e-8; %store E=0 as 1e-8 to be able to take log
      %dk.data{rind}.arefr(recind)  = ...
      %    log10(max(xref_l,abs(dk.data{rind}.xref(recind)*fnrm_e)));  
      
      %error bounds: g11 = g11_i +- g11_e

      dk.data{rind}.nofmn(recind) =line2(1); %number of mn
      dk.data{rind}.lalpha(recind)=line2(2); %number of legendre pols
      dk.data{rind}.g11_e(recind)=line2(3);  %error bound
      dk.data{rind}.g13_e(recind)=line2(4);  %error bound
      dk.data{rind}.g33_e(recind)=line2(5);  %error bound

      %Do not allow positive g11!!
      if dk.data{rind}.g11_i(recind)>0
        if dk.data{rind}.g11_i(recind)-dk.data{rind}.g11_e(recind)<0
          dk.data{rind}.g11_i(recind)=...
              (dk.data{rind}.g11_i(recind)-dk.data{rind}.g11_e(recind))/2;
          dk.data{rind}.g11_e(recind)=abs(dk.data{rind}.g11_i(recind));
        else
          dk.data{rind}.g11_i(recind)=NaN;
        end
      end
      linestr=fgetl(fid); %read beginning of next record
      if length(linestr)>=4
        if strcmp(linestr(2:4),'>3a') || strcmp(linestr(2:4),'c3a')
          %skip the 3a line
          linestr=fgetl(fid);  %read beginning of next record
        end        
      end
      if length(linestr)>=2
        while linestr(2)=='c' %skip the c line
          linestr=fgetl(fid);  %read beginning of next record
        end        
      end
    end
    dk.data{rind}.nrec=recind; %number of cmul,E pairs for this radius
    
    linestr=fgetl(fid); %read first header line for next radius
    if length(linestr)==1
      if linestr==-1 %End of file has been reached
        eof=1;
      end
    end

  end

  dk.nr=rind; %total number of radii

  fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totlen=0;
for rind=1:dk.nr
  
  totlen=totlen+dk.data{rind}.nrec;
  tmp=find(diff(dk.data{rind}.cmul));
  NEH(rind)=tmp(1);                      %number of E values for the runs with the highest cmul
  NEL(rind)=dk.data{rind}.nrec-tmp(end); %number of E values for the runs with the
                                         %lowest cmul
  
  strangecorrfactor=0.1;
  
  g11factor(rind)=-4/3/dk.iota(rind)^2*dk.NormTorCurv(rind)^2;
  invcmultilde=3*dk.data{rind}.EovervB(1:NEH(rind))/dk.iota(rind)^2/dk.r(rind)*dk.R(rind)^2;
  invcmultilde=invcmultilde*strangecorrfactor;
  approxinds= find(dk.data{rind}.cmul(1:NEH(rind)).*invcmultilde <1 & invcmultilde~=0);
  xcmul{rind}=sort(1./invcmultilde(approxinds));
  Nxcmul(rind)=length(xcmul{rind});

  totlen=totlen+(Nxcmul(rind)+1)*NEH(rind)+NEL(rind);
end

extrapolpow1=100; %towards higher collisionalities
extrapolpow1h=50;
extrapolpow2=100; %towards lower collisionalities

r=zeros(totlen,1);
logcmul=zeros(totlen,1);
EovervB=zeros(totlen,1);
g11u=zeros(totlen,1);
g11l=zeros(totlen,1);
g13u=zeros(totlen,1);
g13l=zeros(totlen,1);
g33u=zeros(totlen,1);
g33l=zeros(totlen,1);

thisind=0;
for rind=1:dk.nr
  nE1=NEH(rind);
  nE2=NEL(rind);
  nrec=dk.data{rind}.nrec;
  
  entry_len=nrec+(Nxcmul(rind)+1)*nE1+nE2;
  r(thisind+(1:entry_len)) = dk.r(rind);
  
  logcmul(thisind+(1:nE1))=...
      log(dk.data{rind}.cmul(1:nE1))'+extrapolpow1;
  EovervB(thisind+(1:nE1))=dk.data{rind}.EovervB(1:nE1)';
  for xind=1:Nxcmul(rind)
    logcmul(thisind+xind*nE1+(1:nE1))=log(xcmul{rind}(xind));
  end
  logcmul(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))=log(dk.data{rind}.cmul)';
  logcmul(thisind+(Nxcmul(rind)+1)*nE1+nrec+(1:nE2))=...
      log(dk.data{rind}.cmul(nrec-nE2+1:nrec))'-extrapolpow2;
  
  EovervB(thisind+(1:nE1))=dk.data{rind}.EovervB(1:nE1)';
  for xind=1:Nxcmul(rind)
    EovervB(thisind+xind*nE1+(1:nE1))=dk.data{rind}.EovervB(1:nE1)';
  end
  EovervB(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))=dk.data{rind}.EovervB';
  EovervB(thisind+(Nxcmul(rind)+1)*nE1+nrec+(1:nE2))=dk.data{rind}.EovervB(nrec-nE2+1:nrec)';
  
  %To know how to extrapolate g11 to high collisionalities, we check whether
  % check = 3 vE^* nu^* /(iota eps_t) is < or > 1 for the highest given data.
  invcmultilde=3*dk.data{rind}.EovervB(1:nE1)/dk.iota(rind)^2/dk.r(rind)*dk.R(rind)^2;
  invcmultilde=invcmultilde*strangecorrfactor;
  %cmultilde=dk.iota(rind)^2*dk.r(rind)/dk.R(rind)^2./dk.data{rind}.EovervB(1:nE1)/3;
  %extrapinds= find(dk.data{rind}.cmul(1:nE1)./cmultilde > 1);

  %approxinds= find(dk.data{rind}.cmul(1:nE1).*invcmultilde <1 & invcmultilde~=0);
  %logcmul(thisind+nE1+(approxinds))=-log(invcmultilde(approxinds));

  %g11factor=-4/3/dk.iota(rind)^2*dk.NormTorCurv(rind)^2;
  
  %Extrapolation to high collisionalities
  g11u(thisind+(1:nE1))   =g11factor(rind)*exp(logcmul(thisind+(1:nE1)))...
      ./(1+(exp(logcmul(thisind+(1:nE1))).*invcmultilde').^2);
  g11l(thisind+(1:nE1)) = g11u(thisind+(1:nE1));
  for xind=1:Nxcmul(rind)
    g11u(thisind+xind*nE1+(1:nE1))   =...
        g11factor(rind)*exp(logcmul(thisind+xind*nE1+(1:nE1)))...
        ./(1+(exp(logcmul(thisind+xind*nE1+(1:nE1))).*invcmultilde').^2);
    g11l(thisind+xind*nE1+(1:nE1)) = g11u(thisind+xind*nE1+(1:nE1));
  end
  
  g13u(thisind+(1:(Nxcmul(rind)+1)*nE1))   =zeros((Nxcmul(rind)+1)*nE1,1);
     %0*...
     %  (dk.data{rind}.g13_i(1:nE1)'+dk.data{rind}.g13_e(1:nE1)');
  g13l(thisind+(1:(Nxcmul(rind)+1)*nE1))   =zeros((Nxcmul(rind)+1)*nE1,1);
     %0*...
     %  (dk.data{rind}.g13_i(1:nE1)'-dk.data{rind}.g13_e(1:nE1)');
  
  g33u(thisind+(1:nE1))   =exp(-extrapolpow1)*...
       (dk.data{rind}.g33_i(1:nE1)'+dk.data{rind}.g33_e(1:nE1)');
  g33l(thisind+(1:nE1))   =exp(-extrapolpow1)*...
       (dk.data{rind}.g33_i(1:nE1)'-dk.data{rind}.g33_e(1:nE1)');
  for xind=1:Nxcmul(rind)
    extrapolpow1hX=logcmul(thisind+xind*nE1+(1:nE1))-log(dk.data{rind}.cmul(1:nE1))';
    g33u(thisind+xind*nE1+(1:nE1))   =exp(-extrapolpow1hX).*...
        (dk.data{rind}.g33_i(1:nE1)'+dk.data{rind}.g33_e(1:nE1)');
    g33l(thisind+xind*nE1+(1:nE1))   =exp(-extrapolpow1hX).*...
        (dk.data{rind}.g33_i(1:nE1)'-dk.data{rind}.g33_e(1:nE1)');
  end

  %Extrapolation to low collisionalities
  %assume sqrt(nu) scaling:
  g11u(thisind+(Nxcmul(rind)+1)*nE1+nrec+(2:nE2))   = -ramp( -exp(-extrapolpow2/2)*...
       (dk.data{rind}.g11_i(nrec-nE2+2:nrec)'+dk.data{rind}.g11_e(nrec-nE2+2:nrec)'));
  g11l(thisind+(Nxcmul(rind)+1)*nE1+nrec+(2:nE2))   = -ramp( -exp(-extrapolpow2/2)*...
       (dk.data{rind}.g11_i(nrec-nE2+2:nrec)'-dk.data{rind}.g11_e(nrec-nE2+2:nrec)'));
  %except for the E=0 case, where g11 propto 1/nu
  g11u(thisind+(Nxcmul(rind)+1)*nE1+nrec+1)   =exp(extrapolpow2)*...
       (dk.data{rind}.g11_i(nrec-nE2+1)'+dk.data{rind}.g11_e(nrec-nE2+1)');
  g11l(thisind+(Nxcmul(rind)+1)*nE1+nrec+1)   =exp(extrapolpow2)*...
       (dk.data{rind}.g11_i(nrec-nE2+1)'-dk.data{rind}.g11_e(nrec-nE2+1)');
  
  g13u(thisind+(Nxcmul(rind)+1)*nE1+nrec+(1:nE2))   =...
       (dk.data{rind}.g13_i(nrec-nE2+1:nrec)'+dk.data{rind}.g13_e(nrec-nE2+1:nrec)');
  g13l(thisind+(Nxcmul(rind)+1)*nE1+nrec+(1:nE2))   =...
       (dk.data{rind}.g13_i(nrec-nE2+1:nrec)'-dk.data{rind}.g13_e(nrec-nE2+1:nrec)');

  g33u(thisind+(Nxcmul(rind)+1)*nE1+nrec+(1:nE2))   =exp(extrapolpow2)*...
       (dk.data{rind}.g33_i(nrec-nE2+1:nrec)'+dk.data{rind}.g33_e(nrec-nE2+1:nrec)');
  g33l(thisind+(Nxcmul(rind)+1)*nE1+nrec+(1:nE2))   =exp(extrapolpow2)*...
       (dk.data{rind}.g33_i(nrec-nE2+1:nrec)'-dk.data{rind}.g33_e(nrec-nE2+1:nrec)');

  %rind
  %max(dk.data{rind}.g13_i)
  %max(dk.data{rind}.g13_i'+dk.data{rind}.g13_e')
  %max(dk.data{rind}.g13_i'-dk.data{rind}.g13_e')
  
  g11u(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))   =-ramp(-dk.data{rind}.g11_i'+dk.data{rind}.g11_e');
  g11l(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))   =-ramp(-dk.data{rind}.g11_i'-dk.data{rind}.g11_e');
  g13u(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))   =(dk.data{rind}.g13_i'+dk.data{rind}.g13_e');
  g13l(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))   =(dk.data{rind}.g13_i'-dk.data{rind}.g13_e');
  g33u(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))   =dk.data{rind}.g33_i'+dk.data{rind}.g33_e';
  g33l(thisind+(Nxcmul(rind)+1)*nE1+(1:nrec))   =dk.data{rind}.g33_i'-dk.data{rind}.g33_e';
  
  thisind=thisind+entry_len;  %go to the next entry
end


if not(strcmp(dk.Einterp_type,'no interp'))
  dk.interp.r_factor=dk.r(end)/5;

  if isempty(find(EovervB~=0)) 
    dk.Einterp_type='Er=0';
    disp(['All Er=0 in the file ',filename,' !'])
  end
  
  if strcmp(dk.Einterp_type,'Er=0')
    valids=find(EovervB==0);
    X=[r(valids)/dk.interp.r_factor,logcmul(valids)];
    dk.interp.Flgg11u=TriScatteredInterp(X, log(-g11u(valids)));
    dk.interp.Flgg11l=TriScatteredInterp(X, log(-g11l(valids)));
    dk.interp.Fg13u=TriScatteredInterp(X, g13u(valids));
    dk.interp.Fg13l=TriScatteredInterp(X, g13l(valids));
    dk.interp.Flgg33u=TriScatteredInterp(X, log(-g33u(valids)));
    dk.interp.Flgg33l=TriScatteredInterp(X, log(-g33l(valids)));
  else
    if isempty(find(EovervB==0))
      error('No E=0 entries!')
    end
    dk.interp.EovervBmin=max(EovervB(find(EovervB==0)+1));
    if strcmp(dk.Einterp_type,'asinh')
      Enorm=asinh(EovervB/dk.interp.EovervBmin);
    elseif strcmp(dk.Einterp_type,'linear')
      Enorm=EovervB;
    else
      error('The given dk.Einterp_type is not implemented !')
    end
    X=[r/dk.interp.r_factor,logcmul,Enorm];
    dk.interp.Flgg11u=TriScatteredInterp(X, log(-g11u));
    dk.interp.Flgg11l=TriScatteredInterp(X, log(-g11l));
    dk.interp.Fg13u=TriScatteredInterp(X, g13u);
    dk.interp.Fg13l=TriScatteredInterp(X, g13l);
    dk.interp.Flgg33u=TriScatteredInterp(X, log(-g33u));
    dk.interp.Flgg33l=TriScatteredInterp(X, log(-g33l));
  end  
end


% A help function
function out=ramp(in)
out=in.*(in>0);