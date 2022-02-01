%Run simulation of transcriptional response to a single light pulse
%Runs full diamond IFFL simulation followed by calculation of transcription
%rate etc.

%% Define model parameters

%TFtot = Activator expression level
%Reptot = Repressor expression level
%kdegProt = Fluorescent protein degradation rate = yeast growth rate
%p = Parameter values

%p(1) -> on rate Activator
%p(2) -> off rate Activator
%p(3) -> on rate Repressor
%p(4) -> off rate Repressor
%p(5) -> basal transcription
%p(6) -> max transcription
%p(7) -> MM constant Activator
%p(8) -> Hill coeff Activator
%p(9) -> MM constant Repressor
%p(10) -> Hill coeff Repressor
%p(11) -> mRNA degradation rate
%p(12) -> translation rate / mRNA

TFtot = 1500;
Reptot = 20000;
kdegProt = 0.007;
p = [0.15385,0.02491,0.0060681,0.34393, 0.012389,30.323,1267,3.2682, 1032.41377178639, 3.26765214321773,0.0421160000000000,0.3698];

%% Define experimental parameters and initial conditions

basalRNA = p(5)/p(11);
basalProt = basalRNA * p(12) / kdegProt;
initial = [0 0 basalRNA basalProt];

Imax = 280; % light intensity during pulse
Imin = 0; % light intensity before and after pulse

%% Run response to single pulse

%Response before pulse
I = Imin;
tspan = 0: 0.1 : 10;
[T1,Y1] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt), tspan, initial);	
initial = Y1(end,:);

%Response during pulse
I= Imax;
tspan = 0: 0.1 : 30;
[T2,Y2] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt), tspan, initial);	
initial = Y2(end,:);

%Response after pulse
I = Imin;
tspan = 0: 0.1 : 80;
[T3,Y3] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt), tspan, initial);	
initial = Y3(end,:);

c = [Y1;Y2(2:length(Y2),:);Y3(2:length(Y3),:)]; %concatenation of simulation results

%% Compute transcriptional rate and hill function values
Act = c(:,1);
Rep = c(:,2);    
tr = p(5) + p(6) .* (c(:,1).^ p(8) ./ (c(:,1).^ p(8) + p(7) ^ p(8))) .* (1 ./(1 + (c(:,2)/p(9)).^p(10))) ; % Calculate transcription rate
fA = (c(:,1).^ p(8) ./ (c(:,1).^ p(8) + p(7) ^ p(8))); 
fR = (1 ./(1 + (c(:,2)/p(9)).^p(10)));

%% Plotting
time = (0:(length(tr)-1))./10;

subplot(3,1,1)
plot(time, Act), xlabel('Time (min)'), ylabel('Act(on)')

subplot(3,1,2)
plot(time, Rep), xlabel('Time (min)'), ylabel('Rep(on)')

subplot(3,1,3)
plot(time, tr), xlabel('Time (min)'), ylabel('Transcr. rate')
