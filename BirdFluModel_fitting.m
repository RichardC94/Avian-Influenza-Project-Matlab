function BirdFluModel_fitting
clc
close all

load BirdFluData2016.txt        % imports the bird flu data
tdata = BirdFluData2016(:,1);   % time in years
qdata = BirdFluData2016(:,2);   % number of cumulative cases

load BirdFluData2016_20.txt     % unused real data

tforward = (0:0.01:15)';        % time span for ODE solutions
tmeasure = [1:50:1001]';         % selects points in the solution for 
    ... comparison with data

format long                     % precision

Hpop = 26370;                   % human population of countries affected
    ... in units of 10^5
Bpop = 16540;                   % poultry population of countries affected
    ... in units of 10^6
Ih0 = 0.0005;                   % initial infected human population 
    ... in units of 10^5
betaB = 0.002;                  % initial estimate for transmission rate 
    ... in poultry population
betaH = 5*10^-8;                % initial estimate for transmission rate 
    ... between birds and humans
Ib0 = 5.5;                        % inital estimate for infected 
    ... bird population
        
k = [betaB betaH Ib0];          % initial values of parameters to be fitted


function dx = SIRSI(t,x,k)       % SIRSI equations    

lambdaB = 8270;                 % poultry birth rate
lambda = 26370/65;              % human birth rate
muB = 1/2;                      % poultry natural death rate
nuB = 36.5;                     % poultry death rate due to AI
mu = 1/65;                      % human natural death rate
nu = 36.5;                      % human death rate due to AI
dx = zeros(6,1);                % assigns zeros to dy

dx(1) = lambdaB-k(1)*x(2)*x(1)-muB*x(1);   % change in susceptible 
    ... bird population
dx(2) = k(1)*x(2)*x(1)-(nuB+muB)*x(2);     % change in infected 
    ... bird population
dx(3) = lambda-k(2)*x(2)*x(3)-mu*x(3);     % change in susceptible 
    ... human population
dx(4) = k(2)*x(3)*x(2)-(nu+mu)*x(4);       % change in infected 
    ... human population
dx(5) = k(2)*x(3)*x(2);                    % cumulative cases (no deaths)
dx(6) = nuB*x(2);                          % total poultry deaths due to AI
end

function cases = Inf_hum(k,tdata)    % solves ODEs to give cumulative 
    ... cases values to be compared with real data

[T X] = ode23s(@(t,x)(SIRSI(t,x,k)),tforward,[Bpop-k(3) k(3) Hpop Ih0 Ih0 0]);  
    % solves ODEs with initial conditions
 cases = X(tmeasure(:),5);      % case values corresponding to real data

end

lb = [0 0 0];                   % sets lsqcurvefit lower bound

for i = 1:5

 [k,resnorm,residual] = lsqcurvefit(@Inf_hum,k,tdata,qdata,lb,[],...
 optimset('TolX',10^(-20),'TolFun',10^(-20))); 
    % uses SSE and assigns parameters that give the minimum amount of error
    ... to k

end
 
[T X] = ode23s(@(t,x)(SIRSI(t,x,k)),tforward,[Bpop-k(3) k(3) Hpop Ih0 Ih0 0]);
    % solve ODEs with new parameter values

lambdaB = 8270;                 % poultry birth rate
lambda = 26370/65;              % human birth rate
muB = 1/2;                      % poultry natural death rate
nuB = 36.5;                     % poultry death rate due to AI
Rb = (lambdaB*k(1))/(muB*(muB+nuB))    % display basic reporduction number 
    ... for bird population

betaB = k(1)                    % display parameter values
betaH = k(2)
Ib0 = k(3)

error = resnorm


cases_data1 = BirdFluData2016(:,2).*10^5;              % change units from 10^5 to 1
cases_data2 = BirdFluData2016_20(:,2).*10^5;
cases_calc = X(:,5).*10^5;             % change units from 10^5 to 1
year = 2005+T;                         % change time to year
year_data1 = 2005+BirdFluData2016(:,1);                % change time to year        
year_data2 = 2005+BirdFluData2016_20(:,1);

figure (1)
plot(year_data1,cases_data1,'r.',year_data2,cases_data2,'kx');       % plot real data
hold on
plot(year,cases_calc,'b-');            % plot curve
xlabel('Year');                        % label axis
ylabel('Cumulative Number of Cases');  
legend('Real Data','Unsed Data', 'Equation Solution')

figure (2)
year_res = [2005:0.5:2015];
bar(year_res, residual);
xlabel('Year');                        % label axis
ylabel('Residual'); 
ylim([-10^-3 10^-3])

end