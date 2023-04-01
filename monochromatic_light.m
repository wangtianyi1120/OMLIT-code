%% cc
clear;clc
% declare all data are in nanometer

%% ########### configuration ##########
wave = [390,420,470,555,630,690];

d_tissue = linspace(20,400,39);
% d_tissue = ones(1,length(wave))*60;
% d_coating = linspace(0,400,41);
d_coating = ones(1,length(wave))*0;
d_tape = 50000;
d_sub = 50000;
d_glass = 170000;

n_cell = ones(1,length(wave)) * (1.419+0.0*1i);   % refractive index of cell
n_tissue  = ones(1,length(wave)) * (2.08+1i*1.0);           % refractive index of tissue/background
n_sub = conj(Si_JAW_model(wave/1000));            % refractive index of substrate
n_glass = 1.55;

inc_angle = 0;

para1 = []
para2 = []
para3 = []
para4 = []
Wpara1 = []
Wpara2 = []
Wpara3 = []
FS = 11

%% coating
coating = 'none'
n= [0.21288, 0.15611, 0.098749, 0.072, 0.072, 0.077143]
k= [1.7406, 2.05, 2.5418, 3.296, 3.9078, 4.3775]
n_coating=n+1i*k
% coating = 'Ag'
% n= [0.21288, 0.15611, 0.098749, 0.072, 0.072, 0.077143]
% k= [1.7406, 2.05, 2.5418, 3.296, 3.9078, 4.3775]
% n_coating=n+1i*k
% coating = 'Al'
% n= [0.301, 0.344, 0.424, 0.6025, 0.836, 1.1105]
% k= [3.7, 4.003, 4.491, 5.322, 6.0403, 6.5617]
% n_coating=n+1i*k
% coating = 'C'
% n= [2.6693, 2.6264, 2.6257, 2.7236, 2.8324, 2.9115]
% k= [1.218, 1.2478, 1.3405, 1.4926, 1.5912, 1.6547]
% n_coating=n+1i*k
% coating = 'CNT'
% n= [1.5573, 1.5537, 1.5575, 1.5691, 1.5733, 1.5797]
% k= [0.43088, 0.42203, 0.42066, 0.42679, 0.44286, 0.46679]
% n_coating=n+1i*k
% coating = 'Cr'
% n= [1.965, 2.1228, 2.501, 3.1873, 3.1452, 3.0624]
% k= [2.7906, 2.9728, 3.235, 3.3245, 3.3124, 3.3744]
% n_coating=n+1i*k
% coating = 'Cu'
% n= [1.2685, 1.1923, 1.1372, 0.67635, 0.0925, 0.07252]
% k= [2.0453, 2.1957, 2.4165, 2.4111, 3.4988, 4.1344]
% n_coating=n+1i*k
% coating= 'Pt'
% n= [1.1756, 0.82454, 0.5785, 0.46412, 0.46568, 0.49291]
% k= [2.9896, 3.3108, 4.0306, 5.1701, 6.0954, 6.8078]
% n_coating=n+1i*k


%% tape
tape = 'D50'
n= [1.6339, 1.6232, 1.6094, 1.5943, 1.5864, 1.5821]
k= [1.12e-05, 6.31e-07, 2.39e-07, 4.175e-07, 3.48e-07, 7.08e-08]
n_tape = n+1i*k
% tape = 'PET'    %first nk is not accurate it is for 400nm
% n= [1.6103, 1.6021, 1.5872, 1.5726, 1.5653, 1.5613]
% k= [2.31e-06, 1.76e-06, 1.35e-06, 1.34e-06, 1.23e-06, 8.72e-07]
% n_tape=n+1i*k
% tape = 'PC'
% n= [1.6339, 1.6232, 1.6094, 1.5943, 1.5864, 1.5821]
% k= [1.12e-05, 6.31e-07, 2.39e-07, 4.175e-07, 3.48e-07, 7.08e-08]
% n_tape=n+1i*k

%% ########### Computation ############
% ###### assembly the vectors #######
rs_tissue = zeros(length(d_tissue),length(d_coating));
rp_tissue = zeros(length(d_tissue),length(d_coating));
Rs_tissue = zeros(length(d_tissue),length(d_coating));
Rp_tissue = zeros(length(d_tissue),length(d_coating));
As_tissue = zeros(length(d_tissue),length(d_coating));
Ap_tissue = zeros(length(d_tissue),length(d_coating));
Ts_tissue = zeros(length(d_tissue),length(d_coating));
Tp_tissue = zeros(length(d_tissue),length(d_coating));
rs_cell = zeros(length(d_tissue),length(d_coating));
rp_cell = zeros(length(d_tissue),length(d_coating));
Rs_cell = zeros(length(d_tissue),length(d_coating));
Rp_cell = zeros(length(d_tissue),length(d_coating));
As_cell = zeros(length(d_tissue),length(d_coating));
Ap_cell = zeros(length(d_tissue),length(d_coating));
Ts_cell = zeros(length(d_tissue),length(d_coating));
Tp_cell = zeros(length(d_tissue),length(d_coating));
II = ones(length(d_tissue),length(d_coating));

contrast = zeros(length(d_tissue),length(d_coating));
M = ["contrast","wavelens","tissue",coating]
N = ["95contrast","wavelens","tissue",coating]

m = 0;
% for i = 1 : length(wave)
for i = 3%[1,3,4,5] %1345
    for j = 1 : length(d_tissue)
        for k = 1 : length(d_coating)

            n_vec_tissue = [1,n_tissue(i),n_coating(i),n_tape(i),n_sub(i),1];
            n_vec_cell   = [1,n_cell(i),n_coating(i),n_tape(i),n_sub(i),1];

            % for TE polarization
            polarization = 0;
            % Tissue
            d_vec = [NaN,d_tissue(j),d_coating(k),d_tape,d_sub,NaN];
            [rs_tissue(j,k),~,Rs_tissue(j,k),Ts_tissue(j,k),As_tissue(j,k)] = multilayer_model(wave(i),d_vec,n_vec_tissue,inc_angle,polarization);
            % Cell
            [rs_cell(j,k),~,Rs_cell(j,k),Ts_cell(j,k),As_cell(j,k)] = multilayer_model(wave(i),d_vec,n_vec_cell,inc_angle,polarization);

            % for TM polarization
            polarization = 1;
            % Tissue
            d_vec = [NaN,d_tissue(j),d_coating(k),d_tape,d_sub,NaN];
            [rp_tissue(j,k),~,Rp_tissue(j,k),Tp_tissue(j,k),Ap_tissue(j,k)] = multilayer_model(wave(i),d_vec,n_vec_tissue,inc_angle,polarization);
            % Cell 
            [rp_cell(j,k),~,Rp_cell(j,k),Tp_cell(j,k),Ap_cell(j,k)] = multilayer_model(wave(i),d_vec,n_vec_cell,inc_angle,polarization);

%             Rs_cell(j,k) = II(j,k)-As_cell(j,k)-Ts_cell(j,k)
%             Rp_cell(j,k) = II(j,k)-Ap_cell(j,k)-Tp_cell(j,k)
%             Rs_tissue(j,k) = II(j,k)-As_tissue(j,k)-Ts_tissue(j,k)
%             Rp_tissue(j,k) = II(j,k)-Ap_tissue(j,k)-Tp_tissue(j,k)
            sum_cell   = Rs_cell(j,k) + Rp_cell(j,k);
%             sum_cell   = 2*II-As_cell(j,k)-Ts_cell(j,k)-Ap_cell(j,k)-Tp_cell(j,k);
            sum_tissue = Rs_tissue(j,k) + Rp_tissue(j,k);
%             sum_tissue = 2*II-As_tissue(j,k)-Ts_tissue(j,k)-Ap_tissue(j,k)-Tp_tissue(j,k);
            contrast(j,k) = (sum_cell - sum_tissue) / sum_tissue;
            
            clc
            m = m + 1;
            disp(['Total number is ',num2str(length(d_tissue)*length(d_coating)*length(wave))])
            disp(['Current number is ',num2str(m)])    
        end
    end
    figure

    %%  line
    c = plot(d_tissue,contrast(:,1))
    set(gca,'FontSize',FS)
    xlabel({'tissue','thickness'},'FontSize',FS)
    ylabel('contrast','FontSize',FS)
    ylim([-2 15])
    title([coating,'-coated ',tape,' tape with ',num2str(wave(i)),'nm illumination'],'FontSize',FS)
    y1=max(abs(contrast(:,1)))
    x1 = find(abs(contrast(:,1))==y1)
    hold on
    best0 = max(contrast(:,1))
    x = find(contrast(:,1)==best0)
    y = 1
    plot(d_tissue(x1),contrast(x1,1),'r.','markersize',30)
    legend({[''],['tissue = ',num2str(d_tissue(x1)),'nm' 10 'best contrast = ',num2str(contrast(x1,1))]},'FontSize',FS,'Location','northeast')

%     %%  plot
%     s = surf(d_coating,d_tissue,contrast)
%     alpha(s,'z')
%     set(gca,'FontSize',FS)
%     xlabel({'coating','thickness'},'FontSize',FS)
%     ylabel({'tissue','thickness'},'FontSize',FS)
%     zlabel('contrast','FontSize',FS)
%     zlim([-2 15])
%     title([coating,'-coated ',tape,' tape with ',num2str(wave(i)),'nm illumination'],'FontSize',FS)
%     hold on
%     z1=max(max(abs(contrast)))
%     [x1 y1] = find(abs(contrast)==z1)
%     bc = 0.95*contrast(x1,y1)
% %     [~,y] = min(abs(contrast(x1,:)-bc))
% %     x=x1
%     plot3(d_coating(y1),d_tissue(x1),contrast(x1,y1),'r.','markersize',30)
% %     plot3(d_coating(y),d_tissue(x),contrast(x,y),'g.','markersize',30)
%     legend({[''],['best contrast = ',num2str(contrast(x1,y1)) 10 'tissue thickness = ',num2str(d_tissue(x1)),'nm' 10 'coating thickness = ',num2str(d_coating(y1)),'nm']},'FontSize',FS,'Location','northeast')%,['0.95*best contrast = ',num2str(contrast(x,y)) 10 'coating = ',num2str(d_coating(y)),'nm']},'FontSize',FS,'Location','best')
% %     para1 = [para1,contrast(x,y)]
% %     para2 = [para2,d_coating(y)]
% %     para3 = [para3,d_tissue(x)]
%     M = cat(1,M,[{num2str(contrast(x1,y1))},{num2str(wave(i))},{num2str(d_tissue(x1))},{num2str(d_coating(y1))}])
% %     N = cat(1,N,[{num2str(contrast(x,y))},{num2str(wave(i))},{num2str(d_tissue(x))},{num2str(d_coating(y))}])

    saveas(gcf,['D:\Users\Tyan\OneDrive - USTC\桌面\MATLAB\omlitmodel\simulation_data\agcrcu\',coating,' on ',tape,' for ',num2str(wave(i)),'nm'],'tif')
    end
% [Nsorted,I] = sort(M(:,1))
% Msorted = M(:,I)
% writematrix(M,['C:\Users\Tianyi Wang\Desktop\matlab\omlit\simulation_data\','monolayer','.txt'],'Delimiter','tab','WriteMode','append')
