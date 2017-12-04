function wout=readVMECstruct(filename)

wout=struct();
if isstr(filename)
  info=ncinfo(filename);
  for vi=1:length(info.Variables)
    wout=setfield(wout,info.Variables(vi).Name,...
                       ncread(filename,info.Variables(vi).Name));
  end
else
  error('Unknown input, wout must be the name of the netcdf file!')
end