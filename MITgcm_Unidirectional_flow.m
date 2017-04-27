% Project 4 
% Unidirectional Flow
% MITgcm Two layer flow

%% Read Data
clear;
NT=53;
data=000;
k=150;

for j=1:NT
    
    U=0;
    T=0;
    U_0=0;
    T_0=0;
    
    U = rdmds('U',data, 'n'); % Velocity [m/s]
    T = rdmds('T',data, 'n'); % Temperature [oC]
  
    U_0=U(:,1,1);
    T_0=T(:,1,1);

    for i=2:100;
    
        A=U(:,1,i);  
        B=T(:,1,i);  

        U_0=[U_0,A];  
        T_0=[T_0,B]; 

    end
    
    U_0=U_0';
    T_0=T_0';
    
    eval(sprintf('U%d=U_0',j));
    eval(sprintf('T%d=T_0',j));
    eval(sprintf('Time%d=data*12.4/60/60',j));

    data=data+k;
    
end


RC = rdmds('RC', 'n'); % Z [m]

Z=RC(:,1,1);

for i=2:100;

    D=RC(:,1,i); 
    Z=[Z,D];

end
%%
Depth = rdmds('Depth', 'n'); 

XC = rdmds('XC', 'n'); % X [m]
XG = rdmds('XG', 'n');

RF = rdmds('RF', 'n'); 

XC= XC./1000; %distance in km


%% Removing values beneath seafloor

distance= XC; 
range=-1*Z;


%% Figures

figure(1)
% U contours


for j=1:NT
Variable=eval(sprintf('T%d',j));
k=[12.5 12.5];
[C h]  = contourf(distance,range,Variable,k);
eval(sprintf('x%d=C(1,:)',j));
eval(sprintf('y%d=C(2,:)',j));
 
% Unique x's and their locations
x=eval(sprintf('x%d',j));
y=eval(sprintf('y%d',j));
[ux,~,idx] = unique(x);
% Accumulate 
ymean = accumarray(idx,y,[],@mean);
movymean=movmean(ymean,3);

eval(sprintf('X%d=ux',j));
eval(sprintf('Eta%d=movymean',j));

end

%%
colormap(brewermap(21,'RdBu')) 
v = VideoWriter('\Users\linarotermund\Desktop\TEST.avi')

%    # Unidirectional
%     Fr_0p25_r_0p25_Hm_1
%     Fr_0p75_r_0p25_Hm_1
%     Fr_1_r_0p25_Hm_0p25
%     Fr_0p15_r_0p35_Hm_1
%     Fr_0p55_r_0p35_Hm_1
%     Fr_1_r_0p35_Hm_1
%     Fr_0p15_r_0p55_Hm_1
%     Fr_0p45_r_0p55_Hm_0p75
%     Fr_0p75_r_0p55_Hm_1
%     Fr_0p25_r_0p75_Hm_0p25
%     Fr_0p75_r_0p75_Hm_0p25


open(v)
for j=1:NT
Variable=eval(sprintf('T%d',j));
h=pcolor(distance,range,Variable); %[R h]=contourf(distance,range,Variable,[0:0.2:15]); 
set(h,'EdgeColor','none');
%k=[-2:0.3:2];
%contourf(distance, range, Variable,k)
axis ij;
%ylim([0 100])
title('Temperature Interface')
c=colorbar; c.Label.String='T [^oC]'; caxis([9 15]);
hold on
ylabel('Depth [m]')
xlabel('Distance [km]');
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',20)
Image(j)=getframe(gca)
writeVideo(v, Image(j))
end

%movie(Image)
close(v)
%%
figure(2)

for j=1:3:NT
    hold on
    X=eval(sprintf('X%d',j));
    Y=eval(sprintf('Eta%d',j));
    plot(X,5*Y-j*25,'b-','LineWidth',3)
    axis ij
    hold on
end  
 hold on
 ylabel('Increasing Time')
 xlabel('Distance [km]');
 title('Interface')
 %title('r_u = 0.1, u = 0.2')
 %axis ij
 set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',20);
 set(gca,'YTickLabel',[])
 xlim([40 170])

 %% Composite Froude Number

for j=1:10:NT
    
    I=0;
    Variable=eval(sprintf('T%d',j));
    k=[12.5 12.5];
    figure(3)
    [C h]  = contourf(distance,range,Variable,k);
    
    interface= C(2,2:end);
    interface_x_values=C(1,2:end);

    for ii=1:80
        [M,II]=min(abs(distance(ii)-interface_x_values));
        I=[I,II];
    end

    I=I(2:end);

    interface_80=interface(I);
    interface_x_values_80=interface_x_values(I);

    U_Upper=eval(sprintf('U%d',j));
    U_Lower=eval(sprintf('U%d',j));

    for ii=1:80; 
        w=find(abs(range-interface_80(ii))<2);
        for nn=1:100
            if nn>w(1)          
                U_Upper(nn,ii)=NaN;
            end
            if nn<w(1)          
                U_Lower(nn,ii)=NaN;
            end
        end
    end

    for ii=1:80; 
        w=find(abs(range-Depth(ii))<2);
        for nn=1:100
            if nn>w(1)          
                U_Upper(nn,ii)=NaN;         
                U_Lower(nn,ii)=NaN;
            end
        end
    end

    
    U_layer1 = nanmean(U_Upper);
    U_layer2 = nanmean(U_Lower);

    d1 = max(range)-interface_80;
    d2 = max(range)-d1;

    g_reduced = 2E-3*9.81;

    Fr1_sq = U_layer1.^2./(g_reduced*d1);
    Fr2_sq = U_layer2.^2./(g_reduced*d2);

    G = sqrt(Fr1_sq + Fr2_sq);
    
    figure(4)
    plot(distance,G+j*0.002,'r-','LineWidth',3)
    
    hold on

end


hold on
ylabel('Increasing Time')
xlabel('Distance [km]');
title('Composite Froude Number, G')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]); set(gca,'fontsize',20);
set(gca,'YTickLabel',[])
xlim([50 170])

close(figure(3))
