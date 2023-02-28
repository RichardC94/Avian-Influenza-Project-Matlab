function SIR_fit

clc
close all

load school_data.txt        % load data
format long                 % precision

day = school_data(:,1);     % time coordinates of data
cases = school_data(:,2);   % number of cases each day

tforward = 3:0.01:20;       % time period ODE solved over
tpoints = [1:100:1101]';    % selects the solution data points 
    ... corresponding to real data points given for comparision

gamma = 0.3;                % initial recovery rate estimate
beta = 0.0025;              % initial infection rate estimate
k = [gamma beta];           % parameters
N = 763;                    % total population
I0 = 25;                    % initial infected population
R0 = 0;                     % initial recovered population
S0 = N-(I0+R0);             % initial susceptible population

function dx = SIR_school(t,x,k)    % SIR equations
     
gamma = k(1);               % current estimate for recovery rate
beta = k(2);                % cureent estimate for infection rate
dx = zeros(3,1);            % assigns zeros to dx

dx(1) = x(2) * x(1) * -beta;    % change in susceptible population
dx(2) = x(2) * x(1) * beta - x(2) * gamma; % change in infected population
dx(3) = gamma * x(2);           % change in removed population

end

[T X] = ode23s(@(t,x)(SIR_school(t,x,k)),tforward,[S0 I0 R0]);
                            % solves ODE for plot
Ie = X(:,2);                % infected population estimate

figure(1)
subplot(2,2,1);
plot(day,cases,'r.');       % plot real data
hold on
plot(T,Ie,'b-');            % plot infected population estimate
    ... using initial parameter values 
ylabel('Number of Cases');  % label axis
xlabel('Time in Days');    
legend('Real Data','SIR Model');  
axis([3 20 0 500]);         % axis limits
 
function error = SSE(k)     % error sum of squares

[T X] = ode23s(@(t,x)(SIR_school(t,x,k)),tforward,[S0 I0 R0]); 
                            % solves ODE for SSE

cases_estimate = X(tpoints(:),2); % infected population estimate 

error = sum((cases_estimate - cases).^2); % calculates error sum of squares
    ... between estimated and real data
end

[k,min_error] = fminsearch(@SSE,k);  % adjust parameter values to find 
    ... minimum value for error and assign final paramaters to k

gamma = k(1)                % display final recovery rate
beta = k(2)                 % display final infection rate
R = N*(k(2)/k(1))           % display reproduction number


[T X] = ode23s(@(t,x)(SIR_school(t,x,k)),tforward,[S0 I0 R0]); % solve ODE
    ... with new parameter values

If = X(:,2);                % final infected population estimate

subplot(2,2,2)
plot(day,cases,'r.');       % plot real data
hold on
plot(T,If,'b-');            % plot estimated value of infected population
    ... with final parameter values
ylabel('Number of Cases');  % label axis
xlabel('Time in Days');     
legend('Real Data','SIR Model');    % legend
axis([3 20 0 500]);         % axis limits

subplot(2,2,3)
plot(T,X(:,1),'b-',T,X(:,2),'r',T,X(:,3),'g');   % plots of
    ... susceptible, infected and removed populations against time
xlabel('Days'); ylabel('Number of People');     % label axis
legend('S','I','R');        % legend
axis([3 20 0 800]);         % axis limits

SF = X(:,1)./N;             % susceptible population fraction
IF = X(:,2)./N;             % infected population fraction

subplot(2,2,4)
plot(SF,IF,'r');            % plot phase plane
xlabel('Susceptible Population Fraction');  % label axis
ylabel('Infected Population Fraction');
axis([0 1 0 0.6]);          % axis limits
end