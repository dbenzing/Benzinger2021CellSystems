%Run simulation of transcriptional response to multiple light pulses / PWM
%Runs full diamond IFFL simulation followed by calculation of transcription
%rate. Shown in Supplementary Fig. 2

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

Imax = 280; % light intensity during pulse
Imin = 0; % light intensity before and after pulse

pulsenumber = 5; %number of light pulses simulated
dc = [5,25,50,75]; %Duty-cycles used in simulations
period = 30; %PWM period

basalRNA = p(5)/p(11);
basalProt = basalRNA * p(12) / kdegProt;
res = zeros(length(dc), pulsenumber*period*10+1); %storing results

%% Run pulse responses for different duty-cycles

for i = 1:length(dc)

    initial = [0 0 basalRNA basalProt];
    ontime = dc(i) * period / 100;
    offtime = period - ontime;
    c = initial;
    
    for pi = 1:pulsenumber

        I= Imax;
        tspan = 0: 0.1 : ontime;
        [T1,Y1] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt), tspan, initial);	
        initial = Y1(end,:);

        I = Imin;
        tspan = 0: 0.1 : offtime;
        [T2,Y2] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt), tspan, initial);	
        initial = Y2(end,:);

        c = [c;Y1(2:length(Y1),:);Y2(2:length(Y2),:)];
    end
    
    % Compute transcriptional rate and hill function values  
    res(i,:) = (p(5) + p(6) .* (c(:,1).^ p(8) ./ (c(:,1).^ p(8) + p(7) ^ p(8))) .* (1 ./(1 + (c(:,2)/p(9)).^p(10))))' ; 
end

%% Plotting

time = 0:0.1:(pulsenumber * period);

subplot(4,1,1)
plot(time, res(1,:)), xlabel('Time (min)'), ylabel('Transcr. rate')

subplot(4,1,2)
plot(time, res(2,:)), xlabel('Time (min)'), ylabel('Transcr. rate')

subplot(4,1,3)
plot(time, res(3,:)), xlabel('Time (min)'), ylabel('Transcr. rate')

subplot(4,1,4)
plot(time, res(4,:)), xlabel('Time (min)'), ylabel('Transcr. rate')


