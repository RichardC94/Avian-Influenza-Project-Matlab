function SIR_endemic_plot

clc
close all
format long

N = 50000;                      % total population
tperiod = [0:0.001:5];          % time period
R0 = 0;                         % initial removed population
I0 = 10000;                         % initial infected population
S0 = N-(I0+R0);                 % initial susceptible population

beta = 0.000477;                 % rate of infection
gamma = 26;                     % rate of recovery
mu = 1/2;                       % natural death rate
lambda = N*mu;                  % birth rate

function dx = sir(t, x)         % SIR ODEs

dx = [0; 0; 0];                 % assigns zeroes to dx
beta = 0.000477;                 % rate of infection
gamma = 26;                     % rate of recovery
mu = 1/2;                       % natural death rate
N = 50000;                      % total population
lambda = N*mu;                  % birth rate


dx(1) = lambda - beta * x(1) * x(2) - mu * x(1);            % change in
    ... susceptible population
dx(2) = beta * x(1) * x(2) - gamma * x(2) - mu * x(2);      % change in 
    ... infected population
dx(3) = gamma * x(2) - mu * x(3);   % change in recovered population

end

[t,x] = ode45(@sir, tperiod, [S0 I0 R0]);       % solves ODE

figure (1)
plot(t,x(:,1),'b',t,x(:,2),'r',t,x(:,3),'g','linewidth',1) % plot solutions
xlabel('Time in Years'); ylabel('Population')           % label axis
legend('S','I','R')                             % legend


R = (lambda*beta)/((gamma+mu)*mu)       % basic reproduction number
eqm=[lambda/(mu*R) (mu/beta)*(R-1)]       % endemic equilibrium values (S, I)


for I0 = [1,10000:10000:50000]  % initial infected populations for 
    ... phase portrait

    S0 = N-(I0+R0);             % initial susceptible population

    [t,x] = ode23s(@sir, tperiod, [S0 I0 R0]); % solves ODE

    SF = x(:,1)./N;             % susceptible population fraction
    IF = x(:,2)./N;             % infected population fraction
    T=1-[0:0.1:1];                     % triangle indicating S + I < 1

    figure (2)
    plot(SF,IF,'b',[0:0.1:1],T,'k','linewidth',.5); ylim([0 1]);   % plot 
    ... phase plane
    xlabel('Susceptible Population Fraction'); % label axis
    ylabel('Infected Population Fraction');
    hold on
    plot([0:0.1:1],T,'k')
    plot([0:0.1:1],zeros(1,11),'k')
    plot(zeros(1,11),[0:0.1:1],'k')
end
 
for R0 = 10000:10000:30000      % initial removed populations for phase plane

    I0 = 1;                     % initial infected population    
    S0 = N-(I0+R0);             % initial susceptible population

    [t,x] = ode45(@sir, tperiod, [S0 I0 R0]);   % solves ODE

    SF = x(:,1)./N;             % susceptible population fraction
    IF = x(:,2)./N;             % infected population fraction

    plot(SF,IF,'b','linewidth',.5)  % plot phase plane
    hold on
end
plot(1,0,'kx')
legend('Trajectory','T')        % phase portait legend
plot(eqm(1)./N,eqm(2)./N,'kx','DisplayName','Fixed Point')  % plot fixed 
    ... point
    hold on

end

