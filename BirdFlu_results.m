function BirdFlu_results

clc
close all
format long                     % precision

tforward = (0:0.01:60)';        % time span for ODE solutions
                
Hpop = 26370;                   % human population of countries affected
    ... in units of 10^5
Bpop = 16540;                   % poultry population of countries affected
    ... in units of 10^6
Ih0 = 0.0005;                   % initial infected human population 
    ... in units of 10^5
betaB = 0.002241723385471;      % transmission rate in poultry population
betaH = 0.00000001814169;       % transmission rate between birds and humans
Ib0 = 2.20474587430784;         % infected bird population
        
k = [betaB betaH Ib0];          % fitted parameters
x_lim=2065;

function dx = SISI(t,x,k)       % SISI equations    

lambdaB = 8270;                 % poultry birth rate
lambda = 26370/65;              % human birth rate
muB = 1/2;                      % poultry natural death rate
nuB = 36.5;                     % poultry death rate due to AI
mu = 1/65;                      % human natural death rate
nu = 36.5;                      % human death rate due to AI
dx = zeros(7,1);                % assigns zeros to dy

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
dx(7) = (nuB+muB)*x(2) - muB*x(7);
end
 
[T X] = ode23s(@(t,x)(SISI(t,x,k)),tforward,[Bpop-k(3) k(3) Hpop Ih0 Ih0 0 0]);
    % solve ODEs 

lambdaB = 8270;                 % poultry birth rate
lambda = 26370/65;              % human birth rate
muB = 1/2;                      % poultry natural death rate
nuB = 36.5;                     % poultry death rate due to AI
Rb = (lambdaB*k(1))/(muB*(muB+nuB))    % display basic reporduction number 
    ... for bird population


cases_calc = X(:,5).*10^5;              % change units from 10^5 to 1
year = 2005+T;                          % change time to year

figure (1)                              % human cases curve
plot(year,cases_calc,'b-');             % plot curve
xlabel('Year');                         % label axis
ylabel('Cumulative Number of Cases');  
xlim([2005 x_lim])

figure (2)                              % total poultry deaths due to AI
poultry_deaths = X(:,6).*10^6;          % Change units from 10^6 to 1
plot(year,poultry_deaths,'r');          % plot curve
xlabel('Year');                         % label axis
ylabel('Total Number of Poultry Deaths due to AI');
xlim([2005 x_lim])

figure (3)                              % phase plane plot
plot(X(:,1)*10^6,X(:,2)*10^6,'b');      % change units and plot phase plane
xlabel('Susceptible Bird Population');  % label axis
ylabel('Infective Bird Population');
hold on
eqm = [lambdaB/(muB*Rb) (muB/betaB)*(Rb-1)] % calculate fixed point
plot(eqm(1)*10^6,eqm(2)*10^6,'kx');         % plot fixed point
legend('Trajectory', 'Endemic Equilibrium') % legend

figure (4)                              % infected population against time
plot(year,X(:,2)*10^6,'r')              % plot curve
xlabel('Year')                          % label axis
ylabel('Infected Bird Population')
ylim([0 3e6])                           % y axis limit
xlim([2005 x_lim])
Ib_eqm=(muB/betaB)*(Rb-1)*10^6;
hold on
plot([2005:10:2050],Ib_eqm*ones(1,5),'k--')

figure (5)                              % susceptible population against time
plot(year,X(:,1)*10^6,'b')              % plot curve
xlabel('Year')                          % label axis
ylabel('Susceptible Bird Population')
xlim([2005 x_lim])

figure (6)                              % removed population against time
plot(year,X(:,7)*10^6,'g')              % plot curve
xlabel('Year')                          % label axis
ylabel('Removed Bird Population')
xlim([2005 x_lim])

S_eqm=lambda/(k(2)*eqm(2)+1/65);            % susceptible human equilibrium
I_eqm=(k(2)*S_eqm*eqm(2))/(1/65+36.5);      % infected human equilibrium
rho=1/74;                                   % dimensionless quantity
discriminant=(rho*Rb)^2-4*(rho*Rb-rho);     % eigenvalue discriminant
R_eqm=(nuB/betaB)*(Rb-1);                   % removed bird equilibrium
dI=betaH*S_eqm*eqm(2)                       % human infections per year at equilibrium
dR=nuB*eqm(2)                               % bird deaths per year at equilibrium

figure (7)              % nullcline plot
R=Rb;                   % basic reproduction number
x1=[0:0.1:1.2];         % horizontal y-nullcline x values
y1=zeros(1,13);         % horizontal y-nullcline y values
x2=1/R.*ones(1,15);     % vertical y-nullcline x values
y2=[-0.2:0.1:1.2];      % vertical y-nullcline y values
plot(x1,y1,'b--');      % plot horizontal y-nullcline
hold on
plot(x2,y2,'b--');      % plot vertical y-nullcline
hold on
legend('y-nullcline')   % legend
 
x3=[0.01:0.01:1.2];         % x-nullcline x values
y3=(rho/R).*(1-x3)./x3;     % x-nullcline y values
plot(x3,y3,'r--','DisplayName','x-nullcline')   % x-nullcline
ylim([-.5e-4 1e-4])         % axis limits
xlim([0.994 1.002])         
hold on
 
eqm1=[1,0];                 % disease free equilibrium
plot(eqm1(1),eqm1(2),'kx','DisplayName','Disease-Free Equilibrium')
% plot disease free equilibrium
eqm2=[1/R,rho*(1-1/R)];
% endemic equilibrium
plot(eqm2(1),eqm2(2),'ko','DisplayName','Endemic Equilibrium')
% plot endemic equilibrium
plot(X(:,1)/Bpop,X(:,2)/Bpop,'b','DisplayName','Trajectory');
% plot dimensionless susceptible population against infective 
ylim([-0.5e-4 3e-4])    % axis limits
xlim([0.99 1.004])
xlabel('xb')            % label axis
ylabel('yb')

eigenvalue1=((-rho*R)+((rho*R)^2-4*(rho*R-rho))^0.5)/2
eigenvalue2=((-rho*R)-((rho*R)^2-4*(rho*R-rho))^0.5)/2
end