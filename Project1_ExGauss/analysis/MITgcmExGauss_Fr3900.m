%MITgcm Example MITgcm setup for steady stratified flow over Gaussian Bump

% Snapshots at 5 different times
% Fr = 0.13
clear

%% Read Data

U_t1 = rdmds('U',150, 'n'); %n l a
U_t2 = rdmds('U',300, 'n'); %n l a
U_t3 = rdmds('U',900, 'n'); %n l a
U_t4 = rdmds('U',1800, 'n'); %n l a
U_t5 = rdmds('U',7800, 'n'); %n l a

Depth = rdmds('Depth', 'n'); %n l a

XC = rdmds('XC', 'n'); %n l a %distance in m
XG = rdmds('XG', 'n'); %n l a

XC= XC./1000; %distance in km for T
XG= XG./1000; %distance in km for U

%% Time 
% the time [s] is the file name # multiplied by the time step = 12.4 seconds

Time1=150*12.4/60/60; %Time in hours
Time2=300*12.4/60/60;
Time3=900*12.4/60/60;
Time4=1800*12.4/60/60;
Time5=7800*12.4/60/60;

Str1=sprintf('U [m/s] at time %0.2f hours',Time1)
Str2=sprintf('U [m/s] at time %0.2f hours',Time2)
Str3=sprintf('U [m/s] at time %0.2f hours',Time3)
Str4=sprintf('U [m/s] at time %0.2f hours',Time4)
Str5=sprintf('U [m/s] at time %0.2f hours',Time5)

%% Restructuring U

U_1=U_t1(:,1,1);
U_2=U_t2(:,1,1);
U_3=U_t3(:,1,1);
U_4=U_t4(:,1,1);
U_5=U_t5(:,1,1);

for i=2:25;
    
A=U_t1(:,1,i);
B=U_t2(:,1,i);
C=U_t3(:,1,i);
D=U_t4(:,1,i);
E=U_t5(:,1,i);

U_1=[U_1,A];
U_2=[U_2,B];
U_3=[U_3,C];
U_4=[U_4,D];
U_5=[U_5,E];

end

U_1=U_1.';
U_2=U_2.';
U_3=U_3.';
U_4=U_4.';
U_5=U_5.';

distance= XG; %[1:80];
range=[80:80:2000];

for ii=33:49;%1:length(distance); % Making all values beneath seafloor NaN values
    w=find(abs(Depth(ii)-range')<=80);
    for nn=1:25
        if nn>w
            U_1(nn,ii)=NaN; 
            U_2(nn,ii)=NaN; 
            U_3(nn,ii)=NaN; 
            U_4(nn,ii)=NaN;
            U_5(nn,ii)=NaN;
        end
    end
end

%% Plotting

figure(1)
 
colormap(brewermap(21,'RdBu')) 
 
% Time 1 
subplot(5,1,1);
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_1,k);
set(h,'LineColor','none');
hold on
%xlabel('Distance [km]')
xlim([190 225])
ylabel('Depth [m]')
c=colorbar;
c.Label.String='U [m/s]'
caxis([-3 3])
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
axis ij
title(Str1);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'fontsize',16)

% Time 2
subplot(5,1,2)
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_2,k);
set(h,'LineColor','none');
hold on
%xlabel('Distance [km]')
xlim([190 225])
ylabel('Depth [m]')
shading flat
c=colorbar;
c.Label.String='U [m/s]'
caxis([-3 3])
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
axis ij
title(Str2);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'fontsize',16)

% Time 3
subplot(5,1,3)
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_3,k);
set(h,'LineColor','none');
hold on
%xlabel('Distance [km]')
xlim([190 225])
ylabel('Depth [m]')
shading flat
c=colorbar;
c.Label.String='U [m/s]'
caxis([-3 3])
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
axis ij
title(Str3);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'fontsize',16)

% Time 4
subplot(5,1,4)
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_4,k);
set(h,'LineColor','none');
hold on
%xlabel('Distance [km]')
xlim([190 225])
ylabel('Depth [m]')
shading flat
c=colorbar;
c.Label.String='U [m/s]'
caxis([-3 3])
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
axis ij
title(Str4);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'fontsize',16)

% Time 5
subplot(5,1,5)
k=[-3:0.05:3];
[C h]  = contourf(distance,range,U_5,k);
set(h,'LineColor','none');
hold on
xlabel('Distance [km]')
xlim([190 225])
ylabel('Depth [m]')
shading flat
c=colorbar;
c.Label.String='U [m/s]'
caxis([-3 3])
hold on
p=plot(distance,Depth,'k-','LineWidth',4)
axis ij
title(Str5);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'fontsize',16)

