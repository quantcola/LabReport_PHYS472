dataset=xlsread('data.xlsx');
Bfield=dataset(:,1);    %Tesla
hallvolts=dataset(:,2); %volt
longvolts=dataset(:,3); %volt
resvolts=dataset(:,4);  %volt
resistance=56000;        %ohm
current=resvolts./resistance;   % ampere
realHallVolts=hallvolts+0.6*longvolts;  % the corrected data (eliminate the admixing effect)

HallResistance=realHallVolts./current;
LongResistance=longvolts./current;

HallResistivity=HallResistance;
LongResistivity=LongResistance;

RH=HallResistance./Bfield;

figure(1)
yyaxis left
plot(Bfield,LongResistivity);
hold on
yyaxis right
plot(Bfield,HallResistivity);
legend('Longitudinal Resistivity \rho_{xx}','Hall Resistivity \rho_{xy}')
xlabel("B Field")

figure(2)
plot(Bfield,RH);
legend("Hall Coefficient R_H")
xlabel("B Field (T)")
ylabel("Hall Coefficient R_H")
