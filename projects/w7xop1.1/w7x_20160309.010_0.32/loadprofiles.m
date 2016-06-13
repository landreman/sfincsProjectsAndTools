function [rnorm,r,ne20,ni20,TekeV,TikeV,...
          ErkVm,Zeff,flux21,QMW,comment]...
    =loadprofiles()

filename='w7x_20160309.010_0.32.dat';
fid=fopen(filename);
Ncomment=0;
str=fgetl(fid);
while str(1)=='c' || str(1)=='C'
  Ncomment=Ncomment+1;
  str=fgetl(fid);
end
fclose(fid);
fid=fopen(filename);
for ind=1:Ncomment
  comment{ind}=fgetl(fid);
end
ind=0;
str=fgetl(fid);
while ischar(str)
  ind=ind+1;
  A=sscanf(str,'%f');
  rnorm(ind) =A(1);
  r(ind)     =A(2);
  ne20(ind)  =A(3);
  ni20(ind)  =A(4);
  TekeV(ind) =A(5);
  TikeV(ind) =A(6);
  ErkVm(ind) =A(7);
  Zeff(ind)  =A(8);
  flux21(ind)=A(9); %10^21/s
  QMW(ind)   =A(10);
  
  str=fgetl(fid);
end
fclose(fid);

fig(1);
% Create a plot with 2 y axes using the plotyy function
[ax, h1, h2] = plotyy(rnorm, ne20, rnorm, TekeV, 'plot');
% Add title and x axis label
xlabel('\rho_{tor}');
% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'n_e [10^{20} m^{-3}]');
set(get(ax(2), 'Ylabel'), 'String', 'T [keV]');
axes(ax(2));
hold on 
plot(rnorm,TikeV,'g--')
hold off
legend('T_e','T_i')

fig(2)
plot(rnorm,ErkVm)
ylabel('E_r [kV/m]')
xlabel('\rho_{tor}');

fig(3)
plot(rnorm,flux21*10)
ylabel('\Gamma [10^{20}/s]')
xlabel('\rho_{tor}');

fig(4)
plot(rnorm,QMW)
ylabel('Q [MW]')
xlabel('\rho_{tor}');

