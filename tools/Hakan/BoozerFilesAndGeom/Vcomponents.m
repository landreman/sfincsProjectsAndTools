function [kappa,varpi]=Vcomponents(Booz,vpar)

FSAu        =sum(sum(Booz.u.*Booz.h))/sum(sum(Booz.h));
FSAvparB    =sum(sum(vpar.*Booz.B.*Booz.h))/sum(sum(Booz.h));
FSAvparoverB=sum(sum(vpar./Booz.B.*Booz.h))/sum(sum(Booz.h));


varpi=(FSAvparB/Booz.FSAB2-FSAvparoverB)/FSAu;
kappa=FSAvparB/Booz.FSAB2 - varpi*Booz.G/Booz.FSAB2;

