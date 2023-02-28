function SIR_plot

clc
close all

N = 10000;                  % total population
tperiod = [0:0.01:35];      % time period
R0 = 0;                     % initial removed population
I0 = 1000;                    % initial infected population
S0 = N-(I0+R0);             % initial susceptible population
beta = 0.000045;              % rate of infection
gamma = 0.5;                % rate of recovery    

function dx = SIR(t, x)     % SIR ODEs

dx = [0; 0; 0];           % assigns zeroes to dx
beta = 0.000045;            % rate of infection
gamma = 0.5;              % rate of recovery    

dx(1) = x(2) * x(1) * -beta;      % change in susceptible population
dx(2) = x(2) * beta * x(1) - x(2) * gamma;    % change in infected 
   ...population
dx(3) = x(2) * gamma;             % change in removed population

end

[t,x] = ode45(@SIR, tperiod, [S0 I0 R0]);       % solves ODE

figure (1)
plot(t,x(:,1),'b',t,x(:,2),'r',t,x(:,3),'g'); 
    % plots susceptible, infected and removed populations against time
xlabel('Days'); ylabel('Number of People');     % label axis
legend('S','I','R');

R_number=(beta*N)/gamma

% for phase plane plot a range of initial condition values are needed

for I0 = [10:1000:9010]     % initial infected populations for phase plane

    S0 = N-(I0+R0);             % initial susceptible population

    [t,x] = ode45(@SIR, tperiod, [S0 I0 R0]); % solves ODE

    SF = x(:,1)./N;             % susceptible population fraction
    IF = x(:,2)./N;             % infected population fraction
    T=1-[0:0.1:1];                     % triangle indicating S + I < 1

    figure (2)
    plot(SF,IF,'b',[0:0.1:1],T,'k','linewidth',.5); ylim([0 1]);   % plot 
    ... phase plane
    xlabel('Susceptible Population Fraction'); % label phase plane
    ylabel('Infected Population Fraction');
    hold on
    plot([0:0.1:1],T,'k')
    plot([0:0.1:1],zeros(1,11),'k')
    plot(zeros(1,11),[0:0.1:1],'k')
end

for R0 = 1000:1000:6000     % initial removed populations for phase plane

I0 = 10;                    % initial infected population    
S0 = N-(I0+R0);             % initial susceptible population

[t,x] = ode45(@SIR, tperiod, [S0 I0 R0]);   % solves ODE

SF = x(:,1)./N;             % susceptible population fraction
IF = x(:,2)./N;             % infected population fraction

figure (2);
plot(SF,IF,'b','linewidth',.5)  % plot phase plane
hold on
end

end
