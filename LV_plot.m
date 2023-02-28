function dx = LV_plot(t, x)
clc
close all
function dx = LV(t, x)  % Lotka-Volterra equations
        
dx = [0; 0];    % assign zeros to dx
alpha = 1;      % prey growth rate
beta = .02;     % predation rate
delta = .01;    % relationship between prey population size 
    ... and predator population growth 
gamma = 1;      % predator natural death/migartion rate

dx(1) = - beta * x(1) * x(2) + x(1) * alpha ;  % rate of change of 
    ... prey population
dx(2) = x(2) * x(1) * delta - gamma * x(2); % rate of change of
    ... predator population
                        
end

alpha = 1;                  % prey growth rate
beta = .02;                 % predation rate
delta = .01;                % relationship between prey population size 
    ... and predator population growth 
gamma = 1;                  % predator natural death/migartion rate

y0 = 20;                    % inital predator population
x0 = 20;                    % initial prey population
tperiod = [0:0.01:15];      % time period

[t,x] = ode45(@LV, tperiod, [x0 y0]); % solve ODEs

figure (1)
plot(t,x)               % plot predator and prey populations against time
hold on
xlabel('Time in Months'); ylabel('Population');   % label axis
legend('Prey','Predator');

for y0 = 20:10:40       % initial predator population values 
    ... for phase plane

    [t,x] = ode45(@LV, tperiod, [x0 y0]); % solve ODEs

    figure (2)
    plot(x(:,1),x(:,2),'linewidth',1);    % plot prey population vs predator population 
    ... phase plane
    hold on

end

plot(gamma/delta,alpha/beta,'kx',0,0,'kx')      % plot fixed points
xlabel('Prey Population');      % label phase plane axis
ylabel('Predator Population'); 
legend('y0=20','y0=30','y0=40','Fixed Point');  % legend for initial
    ... predator values

end