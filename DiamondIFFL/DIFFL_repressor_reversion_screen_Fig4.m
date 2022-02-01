% Simulate diamond-IFFL response for different repressor reversion rates

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

TFtot = 2000;
Reptot = 4000;
kdegProt = 0.007;
RepRevRates = 0.1:0.1:1; %Repressor reversion rates used for simulations

%% Define experimental parameters and initial conditions

tspan = [0 360];

basalRNA = p(5)/p(11);
basalProt = basalRNA * p(12) / kdegProt;
initial = [0 0 basalRNA basalProt];

Imax = 280; % PWM light intensity
PWMperiod = 30;
PWMwidth = [0.05:0.05:29.95]; % Width of PWM stimulation, reduce step size for faster simulation

intensities = 0:0.5:280; % AM light intensity, reduce step size for faster simulation

resPWMrev = zeros(length(RepRevRates), length(PWMwidth)); % store PWM results
resAMrev = zeros(length(RepRevRates), length(intensities)); % store AM results

%% Run PWM simulation

r = 1;
for RepRevRate = RepRevRates 
    p = [0.15385,0.02491,0.00606,RepRevRate, 0.012389,30.323,1267,3.2682, 1032.41377178639, 3.26765214321773,0.0421160000000000,0.3698];
    for i = 1:length(PWMwidth)
        resPWMrev(r,i) = PWMsteady(p,initial,TFtot,Reptot, Imax, PWMperiod, PWMwidth(i), kdegProt, tspan);
    end
    r = r + 1;
end

%% Run AM simulation

r = 1;
for RepRevRate = RepRevRates 
    p = [0.15385,0.02491,0.00606,RepRevRate, 0.012389,30.323,1267,3.2682, 1032.41377178639, 3.26765214321773,0.0421160000000000,0.3698];
    ii=1;
    for i = intensities
        [T,Y] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,i,kdegProt), tspan, initial);
        resAMrev(r,ii) = Y(end,end);
        ii = ii +1;
    end
    r = r+1;
end

%% Plotting

subplot(1,2,1)
plot(intensities, resAMrev), xlabel('Light intensity'), ylabel('Protein expression'),title('AM response')

subplot(1,2,2)
plot(PWMwidth, resPWMrev), xlabel('PWM width'), ylabel('Protein expression'),title('PWM response')

