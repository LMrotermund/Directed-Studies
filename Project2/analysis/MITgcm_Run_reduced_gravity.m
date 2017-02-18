% Nonlinear hysdrostatic flow over a gaussian bump
% MITgcm reduced gravity simulation
% Fr= 1.2
% H_m= 0.25

%% Read Data

U_t1 = rdmds('U',150, 'n'); % Velocity [m/s]
U_t2 = rdmds('U',1800, 'n'); 
U_t3 = rdmds('U',7800, 'n'); 

T_t1 = rdmds('T',150, 'n'); % Temperature [^oC]
T_t2 = rdmds('T',1800, 'n'); 
T_t3 = rdmds('T',7800, 'n'); 

Depth = rdmds('Depth', 'n'); 

XC = rdmds('XC', 'n'); % X [m]
XG = rdmds('XG', 'n');

RC = rdmds('RC', 'n'); % Z [m]
RF = rdmds('RF', 'n'); 

XC= XC./1000; %distance in km

%% Time 
% the time [s] is the file name # multiplied by the time step = 12.4 seconds

Time1=150*12.4/60/60; %Time in hours
Time2=1800*12.4/60/60;
Time3=7800*12.4/60/60;

Str1=sprintf('U [m/s] at time %0.2f hours',Time1);
Str2=sprintf('U [m/s] at time %0.2f hours',Time2);
Str3=sprintf('U [m/s] at time %0.2f hours',Time3);
Str4=sprintf('T [^oC] at time %0.2f hours',Time1);
Str5=sprintf('T [^oC] at time %0.2f hours',Time2);
Str6=sprintf('T [^oC] at time %0.2f hours',Time3);

%% Restructuring U, T and Z

U_1=U_t1(:,1,1);
U_2=U_t2(:,1,1);
U_3=U_t3(:,1,1);

T_1=T_t1(:,1,1);
T_2=T_t2(:,1,1);
T_3=T_t3(:,1,1);

Z=RC(:,1,1);

for i=2:50;
    
A=U_t1(:,1,i);  D=U_t2(:,1,i);  E=U_t3(:,1,i);

F=T_t1(:,1,i);  I=T_t2(:,1,i);  J=T_t3(:,1,i);

K=RC(:,1,i);

U_1=[U_1,A];    U_2=[U_2,D];    U_3=[U_3,E];

T_1=[T_1,F];    T_2=[T_2,I];    T_3=[T_3,J];

Z=[Z,K];

end

U_1=U_1.';
U_2=U_2.';
U_3=U_3.';

T_1=T_1.';
T_2=T_2.';
T_3=T_3.';

%% Removing values beneath seafloor

distance= XC; 
range=-1*Z;

for ii=1:80; 
    w=find(abs(Depth(ii)-range')<5);
  
    for nn=1:50
        if nn>w(1)
            U_1(nn,ii)=NaN; 
            U_2(nn,ii)=NaN;
            U_3(nn,ii)=NaN;
            T_1(nn,ii)=NaN; 
            T_2(nn,ii)=NaN;
            T_3(nn,ii)=NaN;
        end
    end
end

%% Figures

figure(1)
 
colormap(brewermap(21,'RdBu')) 
 
% Time 1 
subplot(3,1,1);
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_1,k);
set(h,'LineColor','none');
hold on
ylabel('Depth [m]'); ylim([800 1000])
c=colorbar; c.Label.String='U [m/s]'; caxis([-2 2]);
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
title(Str1);
axis ij
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',16);

% Time 2
subplot(3,1,2)
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_2,k);
set(h,'LineColor','none');
hold on
ylabel('Depth [m]'); ylim([800 1000])
c=colorbar; c.Label.String='U [m/s]'; caxis([-2 2]);
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
title(Str2);
axis ij
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',16);

% Time 3
subplot(3,1,3)
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_3,k);
set(h,'LineColor','none');
hold on
xlabel('Distance [km]'); ylabel('Depth [m]'); ylim([800 1000])
c=colorbar; c.Label.String='U [m/s]'; caxis([-2 2]);
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
title(Str3);
axis ij
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',16);

%%
figure(2)
 
colormap(brewermap(21,'RdBu')) 

% Time 1 
subplot(3,1,1);
k=[0:0.05:10];
[C h]  = contourf(distance,range,T_1,k);
set(h,'LineColor','none');
hold on
ylabel('Depth [m]'); ylim([800 1000])
c=colorbar; c.Label.String='T [^oC]'; caxis([0 10])
hold on
p=plot(distance,Depth-2,'k-','LineWidth',8)
title(Str4);
axis ij
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',16)

% Time 2
subplot(3,1,2)
k=[0:0.05:10];
[C h]  = contourf(distance,range,T_2,k);
set(h,'LineColor','none');
hold on
ylabel('Depth [m]'); ylim([800 1000])
c=colorbar; c.Label.String='T [^oC]'; caxis([0 10])
hold on
p=plot(distance,Depth-2,'k-','LineWidth',8)
title(Str5);
axis ij
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',16)

% Time 3
subplot(3,1,3)
k=[0:0.05:10];
[C h]  = contourf(distance,range,T_3,k);
set(h,'LineColor','none');
hold on
xlabel('Distance [km]'); ylabel('Depth [m]'); ylim([800 1000])
c=colorbar; c.Label.String='T [^oC]'; caxis([0 10])
hold on
p=plot(distance,Depth-2,'k-','LineWidth',8)
title(Str6);
axis ij
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',16)

