%% Import all the data files
FigList = findall(groot, 'Type', 'figure');
for iFig = 1:numel(FigList)
    try
        clf(FigList(iFig));
    catch
        % Nothing to do
    end
end
% (T = 4K) first run
file_field_firstRun_4K=load('data/4K/Field Values -4.20 -1.txt');
file_HallVolts_firstRun_4K=load('data/4K/Hall Volts -4.20 -1.txt');
file_LongVolts_firstRun_4K=load('data/4K/Longitudinal Volts -4.20 -1.txt');
file_ResistorVolts_firstRun_4K=load('data/4K/Resistor Volts -4.20 -1.txt');

% (T = 4K) second run
file_field_secondRun_4K=load('data/4K/Field Values -4.20 -2.txt');
file_HallVolts_secondRun_4K=load('data/4K/Hall Volts -4.20 -2.txt');
file_LongVolts_secondRun_4K=load('data/4K/Longitudinal Volts -4.20 -2.txt');
file_ResistorVolts_secondRun_4K=load('data/4K/Resistor Volts -4.20 -2.txt');

% (T = 4K) third run
file_field_thirdRun_4K=load('data/4K/Field Values -4.20 -3.txt');
file_HallVolts_thirdRun_4K=load('data/4K/Hall Volts -4.20 -3.txt');
file_LongVolts_thirdRun_4K=load('data/4K/Longitudinal Volts -4.20 -3.txt');
file_ResistorVolts_thirdRun_4K=load('data/4K/Resistor Volts -4.20 -3.txt');

% (T = 4K) fourth run
file_field_fourthRun_4K=load('data/4K/Field Values -4.20 -4.txt');
file_HallVolts_fourthRun_4K=load('data/4K/Hall Volts -4.20 -4.txt');
file_LongVolts_fourthRun_4K=load('data/4K/Longitudinal Volts -4.20 -4.txt');
file_ResistorVolts_fourthRun_4K=load('data/4K/Resistor Volts -4.20 -4.txt');

% (T = 1.5K) first run
file_field_firstRun_1_5K=load('data/1.5K/Field Values -4.25 -1.txt');
file_HallVolts_firstRun_1_5K=load('data/1.5K/Hall Volts -4.25 -1.txt');
file_LongVolts_firstRun_1_5K=load('data/1.5K/Longitudinal Volts -4.25 -1.txt');
file_ResistorVolts_firstRun_1_5K=load('data/1.5K/Resistor Volts -4.25 -1.txt');

% The original data was given in '.mat' format, and the format of variable names are
% not uniform, so I output the '.mat' files into '.txt' files before executing this
% analysis program (The txt files were not given at the beginning)
% (T = 1.5K) second run
file_field_secondRun_1_5K=load('data/1.5K/Field Values -4.25 -2.txt');
file_HallVolts_secondRun_1_5K=load('data/1.5K/Hall Volts -4.25 -2.txt');
file_LongVolts_secondRun_1_5K=load('data/1.5K/Longitudinal Volts -4.25 -2.txt');
file_ResistorVolts_secondRun_1_5K=load('data/1.5K/Resistor Volts -4.25 -2.txt');

%% get one-to-one correspondence array of B-field, Hall-Voltage, Longitudinal-Voltage and Resistor-Voltage
% (T = 4K) first run
B_4K_first=file_field_firstRun_4K;
Hall_V_4K_first=(file_HallVolts_firstRun_4K(:,1).^2+file_HallVolts_firstRun_4K(:,2).^2).^0.5;
Long_V_4K_first=(file_LongVolts_firstRun_4K(:,1).^2+file_LongVolts_firstRun_4K(:,2).^2).^0.5;
Res_V_4K_first=(file_ResistorVolts_firstRun_4K(:,1).^2+file_ResistorVolts_firstRun_4K(:,2).^2).^0.5;

% (T = 4K) second run
B_4K_second=file_field_secondRun_4K;
Hall_V_4K_second=(file_HallVolts_secondRun_4K(:,1).^2+file_HallVolts_secondRun_4K(:,2).^2).^0.5;
Long_V_4K_second=(file_LongVolts_secondRun_4K(:,1).^2+file_LongVolts_secondRun_4K(:,2).^2).^0.5;
Res_V_4K_second=(file_ResistorVolts_secondRun_4K(:,1).^2+file_ResistorVolts_secondRun_4K(:,2).^2).^0.5;

% (T = 4K) third run
B_4K_third=file_field_thirdRun_4K;
Hall_V_4K_third=(file_HallVolts_thirdRun_4K(:,1).^2+file_HallVolts_thirdRun_4K(:,2).^2).^0.5;
Long_V_4K_third=(file_LongVolts_thirdRun_4K(:,1).^2+file_LongVolts_thirdRun_4K(:,2).^2).^0.5;
Res_V_4K_third=(file_ResistorVolts_thirdRun_4K(:,1).^2+file_ResistorVolts_thirdRun_4K(:,2).^2).^0.5;

% (T = 4K) fourth run
B_4K_fourth=file_field_fourthRun_4K;
Hall_V_4K_fourth=(file_HallVolts_fourthRun_4K(:,1).^2+file_HallVolts_fourthRun_4K(:,2).^2).^0.5;
Long_V_4K_fourth=(file_LongVolts_fourthRun_4K(:,1).^2+file_LongVolts_fourthRun_4K(:,2).^2).^0.5;
Res_V_4K_fourth=(file_ResistorVolts_fourthRun_4K(:,1).^2+file_ResistorVolts_fourthRun_4K(:,2).^2).^0.5;

% (T = 1.5K) first run
B_1_5K_first=file_field_firstRun_1_5K;
Hall_V_1_5K_first=(file_HallVolts_firstRun_1_5K(:,1).^2+file_HallVolts_firstRun_1_5K(:,2).^2).^0.5;
Long_V_1_5K_first=(file_LongVolts_firstRun_1_5K(:,1).^2+file_LongVolts_firstRun_1_5K(:,2).^2).^0.5;
Res_V_1_5K_first=(file_ResistorVolts_firstRun_1_5K(:,1).^2+file_ResistorVolts_firstRun_1_5K(:,2).^2).^0.5;

% (T = 1.5K) second run
B_1_5K_second=file_field_secondRun_1_5K;
Hall_V_1_5K_second=(file_HallVolts_secondRun_1_5K(:,1).^2+file_HallVolts_secondRun_1_5K(:,2).^2).^0.5;
Long_V_1_5K_second=(file_LongVolts_secondRun_1_5K(:,1).^2+file_LongVolts_secondRun_1_5K(:,2).^2).^0.5;
Res_V_1_5K_second=(file_ResistorVolts_secondRun_1_5K(:,1).^2+file_ResistorVolts_secondRun_1_5K(:,2).^2).^0.5;

%% plot Hall voltages vs B to find when the value of I=(1.6 +- 0.2)e-6 A is measured
%{
figure_1_5K_first=figure(1);
title('\fontsize{16}Hall Voltage and Longitudinal Voltage vs B [\itT]')
xlabel('B [ T ]')
yyaxis left
plot(B_1_5K_first,Hall_V_1_5K_first,'LineWidth',2);
ylabel('\fontsize{15}Hall Voltage [ V ]');
yyaxis right
plot(B_1_5K_first,Long_V_1_5K_first,'LineWidth',2);
ylabel('\fontsize{15}Longitudinal Voltage [ V ]');
saveas(figure_1_5K_first,"figure_1_5K_first.png")
%}

figure_current_candidates=figure(1);
plot(B_1_5K_first,Hall_V_1_5K_first,'LineWidth',2)
hold on
plot(B_1_5K_second,Hall_V_1_5K_second,'LineWidth',2)
hold on
plot(B_4K_first,Hall_V_4K_first,'LineWidth',2)
hold on
plot(B_4K_second,Hall_V_4K_second,'LineWidth',2)
hold on
plot(B_4K_third,Hall_V_4K_third,'LineWidth',2)
hold on
plot(B_4K_fourth,Hall_V_4K_fourth,'LineWidth',2)
legend({'Hall Voltage 1.5K First Run','Hall Voltage 1.5K Second Run','Hall Voltage 4K First Run','Hall Voltage 4K Second Run','Hall Voltage 4K Third Run','Hall Voltage 4K Fourth Run'})
xlabel('B [ T ]')
ylabel('Hall Voltage [ V ]')
saveas(figure_current_candidates,"Fig_current_candidates.png")

%% Estimate the resistance
figure_1_5K_ResV_second=figure(10);
plot(B_1_5K_second,Res_V_1_5K_second,'LineWidth',2);
xlabel('B [ T ]');
ylabel('Resistor Voltage [ V ]');
yline(Res_V_1_5K_second(1),'r--','LineWidth',2)
yline(Res_V_1_5K_second(end),'m-.','LineWidth',2)
legend({'Resistor Voltage at T=4K, Second Run',...
    "Maximum Resistor Voltage = "+num2str(Res_V_1_5K_second(1))+' V',...
    "Lowest Resistor Voltage = "+num2str(Res_V_1_5K_second(end))+' V'})
saveas(figure_1_5K_ResV_second,"figure_1_5K_ResV_second.png");


%% plot the raw data of 4K (Hall voltage and Longitudinal voltage vs B)
figure_4K_first=figure(2);
xlabel('B [ T ]')
yyaxis left
plot(B_4K_first,Hall_V_4K_first,'LineWidth',2);
ylabel('\fontsize{15}Hall Voltage [ V ]');
yyaxis right
plot(B_4K_first,Long_V_4K_first,'LineWidth',2);
ylabel('\fontsize{15}Longitudinal Voltage [ V ]');
saveas(figure_4K_first,"figure_4K_first.png")

figure_4K_second=figure(3);
xlabel('B [ T ]')
yyaxis left
plot(B_4K_second,Hall_V_4K_second,'LineWidth',2);
ylabel('\fontsize{15}Hall Voltage [ V ]');
yyaxis right
plot(B_4K_second,Long_V_4K_second,'LineWidth',2);
ylabel('\fontsize{15}Longitudinal Voltage [ V ]');
saveas(figure_4K_second,"figure_4K_second.png")

figure_4K_third=figure(4);
xlabel('B [ T ]')
yyaxis left
plot(B_4K_third,Hall_V_4K_third,'LineWidth',2);
ylabel('\fontsize{15}Hall Voltage [ V ]');
yyaxis right
plot(B_4K_third,Long_V_4K_third,'LineWidth',2);
ylabel('\fontsize{15}Longitudinal Voltage [ V ]');
saveas(figure_4K_third,"figure_4K_third.png")

figure_4K_fourth=figure(5);
xlabel('B [ T ]')
yyaxis left
plot(B_4K_fourth,Hall_V_4K_fourth,'LineWidth',2);
ylabel('\fontsize{15}Hall Voltage [ V ]');
yyaxis right
plot(B_4K_fourth,Long_V_4K_fourth,'LineWidth',2);
ylabel('\fontsize{15}Longitudinal Voltage [ V ]');
saveas(figure_4K_fourth,"figure_4K_fourth.png")

%% plot how the B field ramps
figure_4K_ramping_first=figure(6);
B_4K_first_index=[1:1:length(B_4K_first)]; %#ok<*NBRAK>
plot(B_4K_first_index,B_4K_first,'LineWidth',2);
xlabel('index (increasing with time)')
ylabel('B [ T ]')
saveas(figure_4K_ramping_first,'figure_4K_ramping_first.png')

figure_4K_ramping_second=figure(7);
B_4K_second_index=[1:1:length(B_4K_second)]; %#ok<*NBRAK>
plot(B_4K_second_index,B_4K_second,'LineWidth',2);
xlabel('index (increasing with time)')
ylabel('B [ T ]')
saveas(figure_4K_ramping_second,'figure_4K_ramping_second.png')

figure_4K_ramping_third=figure(8);
B_4K_third_index=[1:1:length(B_4K_third)]; %#ok<*NBRAK>
plot(B_4K_third_index,B_4K_third,'LineWidth',2);
xlabel('index (increasing with time)')
ylabel('B [ T ]')
saveas(figure_4K_ramping_third,'figure_4K_ramping_third.png')

figure_4K_ramping_fourth=figure(9);
B_4K_fourth_index=[1:1:length(B_4K_fourth)]; %#ok<*NBRAK>
plot(B_4K_fourth_index,B_4K_fourth,'LineWidth',2);
xlabel('index (increasing with time)')
ylabel('B [ T ]')
saveas(figure_4K_ramping_fourth,'figure_4K_ramping_fourth.png')

%% Plot (Hall resistivity and longitudinal resistivity vs B) and (Hall coefficient vs B) at T=4K, first run
resistance=52500;    % ohm
% trim the data array because the measurement is stopped half way
B_4K_first_trimed=B_4K_first(1:416);
Hall_V_4K_first_trimed=Hall_V_4K_first(1:416);
Long_V_4K_first_trimed=Long_V_4K_first(1:416);
Res_V_4K_first_trimed=Res_V_4K_first(1:416);

% minimize the admixing effect
Hall_V_4K_first_trimed=Hall_V_4K_first_trimed+Long_V_4K_first_trimed.*0.8;

% calculate the Hall resistivity and the longitudinal resistivity
current_4K_first=Res_V_4K_first_trimed./resistance;
Hall_resistivity_4K_first=Hall_V_4K_first_trimed./current_4K_first;
Long_resistivity_4K_first=Long_V_4K_first_trimed./current_4K_first;

figure_4K_resistivity_first=figure(11)
yyaxis left
plot(B_4K_first_trimed,Hall_resistivity_4K_first,'LineWidth',2)
ylabel('Hall Resistivity [ \Omega \cdot m ]')
yyaxis right
plot(B_4K_first_trimed,Long_resistivity_4K_first,'LineWidth',2)
% add the plateau index
text(0.4089,1040,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {5}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(0.5115,1040,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)
text(0.6764,1040,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)
text(1.011,1040,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)

ylabel('Longitudinal Resistivity [ \Omega \cdot m ]')
xlabel('B [ T ]')
%add the vertical lines which show the middle of the plateaus
xline(1.011,'m--','LineWidth',2);
xline(0.6764,'m--','LineWidth',2)
xline(0.5115,'m--','LineWidth',2)
xline(0.4089,'m--','LineWidth',2)

saveas(figure_4K_resistivity_first,'figure_4K_resistivity_first.png');

figure_4K_Hall_coefficient_first=figure(12);
Hall_coefficient_4K_first=Hall_resistivity_4K_first./B_4K_first_trimed;
plot(B_4K_first_trimed,Hall_coefficient_4K_first,'LineWidth',1);
xlabel('B [ T ]');
ylabel('Hall Coefficient R_H [ \Omega \cdot m \cdot T^{-1} ]');
saveas(figure_4K_Hall_coefficient_first,'figure_4K_Hall_coefficient_first.png');

%% Plot (Hall resistivity and longitudinal resistivity vs B) and (Hall coefficient vs B) at T=1.5K, first run
% trim the data array because the measurement is stopped half way
B_1_5K_first_trimed=B_1_5K_first(1:663);
Hall_V_1_5K_first_trimed=Hall_V_1_5K_first(1:663);
Long_V_1_5K_first_trimed=Long_V_1_5K_first(1:663);
Res_V_1_5K_first_trimed=Res_V_1_5K_first(1:663);

% minimize the admixing effect
Hall_V_1_5K_first_trimed=Hall_V_1_5K_first_trimed+Long_V_1_5K_first_trimed.*0.8;

% calculate the Hall resistivity and the longitudinal resistivity
current_1_5K_first=Res_V_1_5K_first_trimed./resistance;
Hall_resistivity_1_5K_first=Hall_V_1_5K_first_trimed./current_1_5K_first;
Long_resistivity_1_5K_first=Long_V_1_5K_first_trimed./current_1_5K_first;

figure_1_5K_resistivity_first=figure(13)
yyaxis left
plot(B_1_5K_first_trimed,Hall_resistivity_1_5K_first,'LineWidth',2)
ylabel('Hall Resistivity [ \Omega \cdot m ]')
yyaxis right
plot(B_1_5K_first_trimed,Long_resistivity_1_5K_first,'LineWidth',2)
% add the plateau index
text(2.016,2050,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(1.345,2050,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(1.010,2050,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(0.7973,2050,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {5}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)
text(0.6751,2050,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {6}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)

ylabel('Longitudinal Resistivity [ \Omega \cdot m ]')
xlabel('B [ T ]')
%add the vertical lines which show the middle of the plateaus
xline(2.016,'m--','LineWidth',2);
xline(1.345,'m--','LineWidth',2);
xline(1.010,'m--','LineWidth',2);
xline(0.7973,'m--','LineWidth',2)
xline(0.6764,'m--','LineWidth',2)


saveas(figure_1_5K_resistivity_first,'figure_1_5K_resistivity_first.png');

figure_1_5K_Hall_coefficient_first=figure(14);
Hall_coefficient_1_5K_first=Hall_resistivity_1_5K_first./B_1_5K_first_trimed;
plot(B_1_5K_first_trimed,Hall_coefficient_1_5K_first,'LineWidth',1);
xlabel('B [ T ]');
ylabel('Hall Coefficient R_H [ \Omega \cdot m \cdot T^{-1} ]');
saveas(figure_1_5K_Hall_coefficient_first,'figure_1_5K_Hall_coefficient_first.png');


%% Plot (Hall resistivity and longitudinal resistivity vs B) and (Hall coefficient vs B) at T=1.5K, second run
B_1_5K_second_trimed=B_1_5K_second(1:end);
Hall_V_1_5K_second_trimed=Hall_V_1_5K_second(1:end);
Long_V_1_5K_second_trimed=Long_V_1_5K_second(1:end);
Res_V_1_5K_second_trimed=Res_V_1_5K_second(1:end);

% minimize the admixing effect
Hall_V_1_5K_second_trimed=Hall_V_1_5K_second_trimed+Long_V_1_5K_second_trimed.*0.8;

% calculate the Hall resistivity and the longitudinal resistivity
current_1_5K_second=Res_V_1_5K_second_trimed./resistance;
Hall_resistivity_1_5K_second=Hall_V_1_5K_second_trimed./current_1_5K_second;
Long_resistivity_1_5K_second=Long_V_1_5K_second_trimed./current_1_5K_second;

figure_1_5K_resistivity_second=figure(15)
yyaxis left
plot(B_1_5K_second_trimed,Hall_resistivity_1_5K_second,'LineWidth',2)
ylabel('Hall Resistivity [ \Omega \cdot m ]')
yyaxis right
plot(B_1_5K_second_trimed,Long_resistivity_1_5K_second,'LineWidth',2)
% add the plateau index
text(4.019,2570,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(2.676,2570,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(2.013,2570,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16 )
text(1.607,2570,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {5}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)
text(1.345,2570,'\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {6}}}','Color','m','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter', 'latex','fontsize', 16)

ylabel('Longitudinal Resistivity [ \Omega \cdot m ]')
xlabel('B [ T ]')
%add the vertical lines which show the middle of the plateaus
xline(4.019,'m--','LineWidth',2);
xline(2.676,'m--','LineWidth',2);
xline(2.013,'m--','LineWidth',2);
xline(1.607,'m--','LineWidth',2)
xline(1.345,'m--','LineWidth',2)


saveas(figure_1_5K_resistivity_second,'figure_1_5K_resistivity_second.png');

figure_1_5K_Hall_coefficient_second=figure(16);
Hall_coefficient_1_5K_second=Hall_resistivity_1_5K_second./B_1_5K_second_trimed;
plot(B_1_5K_second_trimed,Hall_coefficient_1_5K_second,'LineWidth',1);
xlabel('B [ T ]');
ylabel('Hall Coefficient R_H [ \Omega \cdot m \cdot T^{-1} ]');
saveas(figure_1_5K_Hall_coefficient_second,'figure_1_5K_Hall_coefficient_second.png');