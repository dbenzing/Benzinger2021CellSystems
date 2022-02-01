%Run simulations of the demultiplexer in response to AM and PWM

%% Define model parameters

%p are the parameters of the diamond-IFFL

%p(1) -> on rate Act
%p(2) -> off rate Act
%p(3) -> on rate Rep
%p(4) -> off rate Rep
%p(5) -> basal transcription
%p(6) -> max transcription
%p(7) -> MM constant Act
%p(8) -> Hill coeff Act
%p(9) -> MM constant Rep
%p(10) -> Hill coeff Rep
%p(11) -> mRNA degradation rate
%p(12) -> translation rate / mRNA
%kdegProt -> Fluorescent protein degradation rate

%p2 are the parameters of the second gene expression system

%p2(1) -> on rate Act
%p2(2) -> off rate Act
%p2(3) -> basal transcription
%p2(4) -> max transcription
%p2(5) -> MM constant Act
%p2(6) -> Hill coeff Act
%p(7) -> mRNA degradation rate
%p(8) -> translation rate / mRNA

TFtot = 2000;
TFtot2 = 2000;
Reptot = 10000;
p = [0.15385,0.02491,0.0039878,0.31324, 0.012389,30.323,1267,3.2682, 1032.41377178639, 3.26765214321773,0.0421160000000000,0.3698];
p2 = [0.0060681,0.277, 0.012389,15,1000,5.2682,0.0421160000000000,0.3698];
kdegProt = 0.007;

%% Define experiment and initial conditions

tspan = [0 360]; % timespan of simulated experiments

Imax = 210; % maximal light intensities used in experiments
Imin = 0; % minimal light intensities used in experiments

basalRNA1 = p(5)/p(11);
basalProt1 = basalRNA1 * p(12) / kdegProt;
basalRNA2 = p2(3)/p2(8);
basalProt2 = basalRNA2 * p2(8) / kdegProt;
initial = [0 0 basalRNA1 basalProt1 0 basalRNA2 basalProt2]; % initial conditions for simulation


%% Simulate AM dose response

resAM = zeros(2,length(0:Imax));

for i = 0:Imax
    [T1,Y1] = ode23s(@(t,y) detDemulti(t,y,p,p2,TFtot,TFtot2,Reptot,i,kdegProt), tspan, initial);
    resAM(1,i+1) = Y1(end,4); %Diamond-IFFL response
    resAM(2,i+1) = Y1(end,7); %Second system response
end


%% Simulate PWM with 30 min period and high light intensity

period = 30;
width=[0.05:0.05:29.95];

resPWM = zeros(2,length(width));

for i = 1: length(width)
    [o1, o2] = MultiPWMsteady(p,p2,TFtot,TFtot2,Reptot,Imax, period, width(i),kdegProt,tspan,initial);
    resPWM(1,i) = o1;
    resPWM(2,i) = o2;
end

%% Simulate PWM with 30 min period and low light intensity

Ilow = 3.5;
period = 30;
width=[0.05:0.05:29.95];

resPWMlow = zeros(2,length(width));

for i = 1: length(width)
    [o1, o2] = MultiPWMsteady(p,p2,TFtot,TFtot2,Reptot,Ilow, period, width(i),kdegProt,tspan,initial);
    resPWMlow(1,i) = o1;
    resPWMlow(2,i) = o2;
end

%% Plotting

subplot(1,3,1)
plot(0:Imax, resAM), xlabel('Light intensity'), ylabel('Protein expression'), title('AM')

subplot(1,3,2)
plot(width, resPWM), xlabel('Light intensity'), ylabel('Protein expression'), title('PWMhigh')

subplot(1,3,3)
plot(width, resPWMlow), xlabel('Light intensity'), ylabel('Protein expression'), title('PWMlow')