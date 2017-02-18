% Directed Studies Project #2

% Recreating Figure 2.6 (Baines) using the method of characteristics

% 1d nonlinear hydrostatic flow over an obstacle centred at 0 and
% approaching a height of 0 of x= +/- 3

%% Setup

clear

x3=[-20:0.4:20];        % x-position

dt=0.2;                 % delta T
t_o=0;                  % initial time
t=[t_o:dt:30];          % time array

NT= length(t);          % number of time steps
NX=length(x3);          % number of x steps

g=9.81;                 % gravity

Fr=0.8; % 1.2; %        % Froude number
u_o=1;                  % initial velocity
d_o=((u_o/Fr).^2)/g;    % initial depth

xchar= [-6:0.3:6];      % x- characteristics

H=0.25*d_o;             % height of obstacle
h=exp(((-x3.^2)./3))*H; % obstacle
h_x=gradient(h);        % gradient of obstacle dh/dx
%h_x=-2*a*exp(-x3.^2).*x3; 

method='linear';        % for interpolation

%% Initializing

Ui = zeros(NT,NX);
Ui(1,:) = u_o;

D = zeros(NT,NX);
D(1,:) = d_o;

char= zeros(NT,length(xchar));
char(1,:)=xchar;

%% Calculating characteristics for all x and t

for ii=2:NT;

    for kk=1:NX;
        
        xx=x3(kk);
        
        u3=Ui(ii-1,kk);
        d3=D(ii-1,kk);

        x1 = xx - (u3 + sqrt(g*d3))*dt;
        x2 = xx - (u3 - sqrt(g*d3))*dt;
        
        d1 = interp1(x3, D(ii-1,:),x1,method,'extrap');
        d2 = interp1(x3, D(ii-1,:),x2,method,'extrap');
        
        u1 = interp1(x3, Ui(ii-1,:),x1,method,'extrap');
        u2 = interp1(x3, Ui(ii-1,:),x2,method,'extrap');
        
        zeta = u1 + 2*sqrt(g*d1);
        eta = u2 - 2*sqrt(g*d2);
            
        dhdx1 = interp1(x3, h_x,x1,method,'extrap');
        dhdx2 = interp1(x3, h_x,x2,method,'extrap');
            
        Ui(ii, kk) = (-g*(dhdx1 + dhdx2)*dt + zeta + eta)/2;
        Ci= ( -g*(dhdx1 - dhdx2)*dt +zeta - eta)/4;
        D(ii, kk) = (Ci.^2)./g;
                            
    end

    q = char(ii-1 , :);
    uq=interp1(x3, Ui(ii, :), q,method,'extrap');
    dq = interp1(x3, D(ii, :), q,method,'extrap');
    char(ii, :) = (uq - sqrt(g*dq))*dt + q;

end

%% Figures

figure(1)
subplot(1,2,1)
plot(x3, h/d_o, 'k','LineWidth',2)
xlabel('X')
ylabel('Height (h/d_o)')
set(gca,'FontSize',24)
xlim([-10 10])

subplot(1,2,2)
plot(x3, h_x, 'k','LineWidth',2)
xlabel('X')
ylabel('dh/dx')
set(gca,'FontSize',24)
xlim([-10 10])

figure(2)

subplot(1,2,1)
plot(char, t, 'k','LineWidth',2)
xlabel('X')
ylabel('Time')
set(gca,'FontSize',24)
title('Froude Number= 0.8')
ylim([0 30])
xlim([-10 10])

subplot(1,2,2)
hold on
nn=-1;
for jj = 1:10:100
 plot(x3,D(jj,:)/d_o+nn,'LineWidth',2)
 nn=nn+0.1;
end
xlabel('X')
ylabel('Increasing Time')
title('Interface Height')
xlim([-10 10])
set(gca,'FontSize',24)
set(gca,'YTick',[]);

