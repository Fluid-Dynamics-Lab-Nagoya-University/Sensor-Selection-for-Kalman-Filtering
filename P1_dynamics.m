b%% Main program
%% //////////////////////////////////////////////////////////////
% Comments:
% 	Collaborator: Kumi Nakai, Taku Nonomura 2020
% 	Last modified: 2021/12/24 K.Nakai
% Nomenclature:
% - Scalars
%   n : Number of degrees of freedom of spatial POD modes (state dimension)
%   p : Number of sensors
%   r : Number of rank for truncated POD
%   s : Number of components
% - Matrices
% 	X : Supervising data matrix
% 	Y : Observation matrix
% 	H : Sparse sensor location matrix
% 	U : Spatial POD modes matrix
% 	C : Measurement matrix
% 	Z : POD mode amplitude matrix
%% ==============================================================

close all; clearvars;
tStart=tic;

%% Selection of Problems =======================================
  num_problem=1; % //Random dynamical system//
% num_problem=2; % //Linear system of three dynamical modes//
% num_problem=3; % //PIV// (code for NOAA-SST has not yet developed) 

%% Input parameters ============================================
pmin = 1;
pinc = 1;
pmax = 10;
ps   = 10;%pmin:pinc:pmax;
num_ave = 3;
if num_problem == 1
    num_ave = 3;
    r = 5;
    n = 1000;
    sigma_s2 = 1e-1; %system noise
    sigma_o2 = 1e0; %observation noise
elseif num_problem == 2
    num_ave = 3;
    r = 6;
    n = 1000;
    f = [1.0 2.5 5.5];
    g = [-0.001 -0.001 -0.03];
    sigma_s2 = 1e-1;%system noise
    sigma_o2 = 1e0;%observation noise
elseif num_problem == 3
    r = 10;
    s = 2;
    test_date = 20210518; run_num = 148; AoA = 18;
end

%% Preparation of output directories ===========================
srcdir  = pwd;
if (num_problem == 1) || (num_problem == 2)
    workdir = ('../work');
    workdir2  = [workdir,'/fig_loop'];
    mkdir(workdir); mkdir(workdir2);
elseif num_problem == 3
    datadir = '../data';
    workdir = '../work_sensor';
    workdir2 = [workdir,'/',num2str(AoA),'deg_s',num2str(s)];
    workdir3 = [workdir2,'/r',num2str(r)];
    if not(exist(workdir,'dir'))
        mkdir(workdir); disp(['New directory "',workdir,'" was made']);
    end
    if not(exist(workdir2,'dir'))
        mkdir(workdir2); disp(['New directory "',workdir2,'" was made']);
    end
    if not(exist(workdir3,'dir'))
        mkdir(workdir3); disp(['New directory "',workdir3,'" was made']);
    end
end

%% Random system problem =======================================
if (num_problem == 1) || (num_problem == 2)
    rng(3,'twister');
    for w=1:1:num_ave
        text = [ 'Average loop: ', num2str(w),'/', num2str(num_ave) ];
        disp(text);
        %% Generation dynamical system ==========================
        if num_problem == 1
            % /// Random dynamical system ///            
            B_sys = 0;
            cond_AAI = 1e10;
            CNT2=0;
            while (B_sys == 0) || (cond_AAI > 1e3)
                CNT2 = CNT2+1;
                sysorg = drss(r,n); % generate random discrete test model
                % (r-th order model with one input and n outputs)
                B_sys = isstable(sysorg); % Boolean value of stability
                [Aorg,Borg,Corg,Dorg] = ssdata(sysorg);
                AAI = kron(Aorg',Aorg')-eye(r^2);
                cond_AAI = cond(AAI);
                text = [ 'CNT2=', num2str(CNT2), ', B_sys=', num2str(B_sys),...
                    ', cond(AAI)=', num2str(cond_AAI) ];
                disp(text);
            end
        else
            % /// Linear system of three dynamical modes ///
            % specify characteristic frequencies and growth/decay rates
            % associated with continuous-time dynamics
            assert(length(f)==length(g))
            % construct low-rank continuous-time operator (rank=k)
            k = 2*length(f);  % (to perform DMD/TDMD with correct rank, set r=k)
            A1 = [];
            for ii = 1:length(f)
                A2 = [[g(ii) 2*pi*f(ii); -2*pi*f(ii) g(ii)]];
                A1 = [A1 A2];
            end
            Alowrank = [];
            for ii = 1:length(f)
                Alowrank = blkdiag(Alowrank,A1(:,(ii-1)*2+1:2*ii));
            end
            % size(Alowrank)
            dt = 0.01; % time step size
            ii = 1;
            Aorg = expm(dt*ii*Alowrank);
            AAI = kron(Aorg',Aorg')-eye(r^2);
            sysorg = drss(r,n);
            [~,Borg,Corg,Dorg] = ssdata(sysorg);
        end
        
        %% Sensor selection =====================================
        CNT=0;
        for p = ps
            CNT = CNT+1;
            text = [ num2str(p),' sensor selection started --->' ];
            disp(text);
                        
            %% Dynamics w/ noise---------------------------------
            % Minimization of steady state error covariance matrix (Kalman filter)
            [time_RAND(CNT,w+1), H_RAND, sensors_RAND ] = F_sensor_random(n,p);
            trWo_RAND   (CNT,w+1) = F_calc_trWo (  H_RAND,Aorg,     Corg);
            trWo_RAND_2 (CNT,w+1) = F_calc_trWo2(p,H_RAND,Aorg,     Corg);
            trWo_RAND_3 (CNT,w+1) = F_calc_trWo3(  H_RAND,Aorg,     Corg);
            trP_RAND    (CNT,w+1) = F_calc_trP  (p,H_RAND,Aorg,Borg,Corg,sigma_s2,sigma_o2);
            trCC_RAND   (CNT,w+1) = F_calc_tr   (p,H_RAND,          Corg);   

            % Weight descrete CVX
            [time_GSKF(CNT,w+1), H_GSKF, sensors_GSKF ] = F_sensor_GSDP_KF(Aorg,Corg,sigma_s2*eye(size(Aorg)),sigma_o2,p);
            trWo_GSKF   (CNT,w+1) = F_calc_trWo (  H_GSKF,Aorg,     Corg);
            trWo_GSKF_2 (CNT,w+1) = F_calc_trWo2(p,H_GSKF,Aorg,     Corg);
            trWo_GSKF_3 (CNT,w+1) = F_calc_trWo3(  H_GSKF,Aorg,     Corg);
            trP_GSKF    (CNT,w+1) = F_calc_trP  (p,H_GSKF,Aorg,Borg,Corg,sigma_s2,sigma_o2);
            trCC_GSKF   (CNT,w+1) = F_calc_tr   (p,H_GSKF,          Corg);   

            % Weight descrete CVX
            [time_WSKF(CNT,w+1), H_WSKF, sensors_WSKF ] = F_sensor_WSDP_KF(Aorg,Corg,sigma_s2*eye(size(Aorg)),sigma_o2,p);
            trWo_WSKF   (CNT,w+1) = F_calc_trWo (  H_WSKF,Aorg,     Corg);
            trWo_WSKF_2 (CNT,w+1) = F_calc_trWo2(p,H_WSKF,Aorg,     Corg);
            trWo_WSKF_3 (CNT,w+1) = F_calc_trWo3(  H_WSKF,Aorg,     Corg);
            trP_WSKF    (CNT,w+1) = F_calc_trP  (p,H_WSKF,Aorg,Borg,Corg,sigma_s2,sigma_o2);
            trCC_WSKF   (CNT,w+1) = F_calc_tr   (p,H_WSKF,          Corg);                        
            
            % Weight Random Convex
            [time_WRCKF(CNT,w+1), H_WRCKF, sensors_WRCKF ] = F_sensor_CRconvex_KF(Aorg,Corg,sigma_s2*eye(size(Aorg)),sigma_o2,p,1000,n/10);
            trWo_WRCKF   (CNT,w+1) = F_calc_trWo (  H_WRCKF,Aorg,     Corg);
            trWo_WRCKF_2 (CNT,w+1) = F_calc_trWo2(p,H_WRCKF,Aorg,     Corg);
            trWo_WRCKF_3 (CNT,w+1) = F_calc_trWo3(  H_WRCKF,Aorg,     Corg);
            trP_WRCKF    (CNT,w+1) = F_calc_trP  (p,H_WRCKF,Aorg,Borg,Corg,sigma_s2,sigma_o2);
            trCC_WRCKF   (CNT,w+1) = F_calc_tr   (p,H_WRCKF,          Corg);            

            % Weight Convex
            [time_WCKF(CNT,w+1), H_WCKF, sensors_WCKF ] = F_sensor_convex_KF(Aorg,Corg,sigma_s2*eye(size(Aorg)),sigma_o2,p,1000);
            trWo_WCKF   (CNT,w+1) = F_calc_trWo (  H_WCKF,Aorg,     Corg);
            trWo_WCKF_2 (CNT,w+1) = F_calc_trWo2(p,H_WCKF,Aorg,     Corg);
            trWo_WCKF_3 (CNT,w+1) = F_calc_trWo3(  H_WCKF,Aorg,     Corg);
            trP_WCKF    (CNT,w+1) = F_calc_trP  (p,H_WCKF,Aorg,Borg,Corg,sigma_s2,sigma_o2);
            trCC_WCKF   (CNT,w+1) = F_calc_tr   (p,H_WCKF,          Corg);            
            
            % Minimization of steady state error covariance matrix (Kalman filter)
            [time_GKF(CNT,w+1), H_GKF, sensors_GKF, ranki(CNT,:)] = F_sensor_GKF(Aorg,Borg,Corg,p,sigma_s2,sigma_o2);
            trWo_GKF   (CNT,w+1) = F_calc_trWo (  H_GKF,Aorg,     Corg);
            trWo_GKF_2 (CNT,w+1) = F_calc_trWo2(p,H_GKF,Aorg,     Corg);
            trWo_GKF_3 (CNT,w+1) = F_calc_trWo3(  H_GKF,Aorg,     Corg);
            trP_GKF    (CNT,w+1) = F_calc_trP  (p,H_GKF,Aorg,Borg,Corg,sigma_s2,sigma_o2);
            trCC_GKF   (CNT,w+1) = F_calc_tr   (p,H_GKF,          Corg);
            
%             %% Dynamics w/o noise -------------------------------
%             % Maximization of obserbility Gramian
%             [time_Gram(CNT,w+1), H_Gram, sensors_Gram] = F_sensor_Gram(Aorg,Corg,p);
%             detWo_Gram   (CNT,w+1) = F_calc_detWo (H_Gram,Aorg,Corg);
%             detWo_Gram_2 (CNT,w+1) = F_calc_detWo2(p,H_Gram,Aorg,Corg);
%             detWo_Gram_3 (CNT,w+1) = F_calc_detWo3(H_Gram,Aorg,Corg);
%             detP_Gram (CNT,w+1)    = F_calc_detP  (p,H_Gram,Aorg,Borg,Corg,sigma_s2,sigma_o2);
%             detCC_Gram (CNT,w+1)   = F_calc_det (p,H_Gram,Corg);
%             
%             %% S2S w/ noise--------------------------------------
%             % Maximization of determinant of (CRC)
%             [time_S2SwR(CNT,w+1), H_S2SwR, sensors_S2SwR] = F_sensor_DGwR(Corg,p,sigma_o2);
%             detWo_S2SwR   (CNT,w+1) = F_calc_detWo (H_S2SwR,Aorg,Corg);
%             detWo_S2SwR_2 (CNT,w+1) = F_calc_detWo2(p,H_S2SwR,Aorg,Corg);
%             detWo_S2SwR_3 (CNT,w+1) = F_calc_detWo3(H_S2SwR,Aorg,Corg);
%             detP_S2SwR (CNT,w+1)    = F_calc_detP  (p,H_S2SwR,Aorg,Borg,Corg,sigma_s2,sigma_o2);
%             detCC_S2SwR (CNT,w+1)   = F_calc_det (p,H_S2SwR,Corg);
%             
%             %% S2S w/o noise-------------------------------------
%             % Maximization of determinant of Fisher information matrix
%             [time_S2S(CNT,w+1), H_S2S, sensors_S2S] = F_sensor_QD(Corg,p);
%             detWo_S2S   (CNT,w+1) = F_calc_detWo (H_S2S,Aorg,Corg);
%             detWo_S2S_2 (CNT,w+1) = F_calc_detWo2(p,H_S2S,Aorg,Corg);
%             detWo_S2S_3 (CNT,w+1) = F_calc_detWo3(H_S2S,Aorg,Corg);
%             detP_S2S (CNT,w+1)    = F_calc_detP  (p,H_S2S,Aorg,Borg,Corg,sigma_s2,sigma_o2);
%             detCC_S2S (CNT,w+1)   = F_calc_det (p,H_S2S,Corg);
            
            text = [ '---> ', num2str(p), ' sensor selection finished!' ];
            % disp(text);
        end
        
        %% Temporary plot ======================================
        cd(workdir2)
        newcolors = {'red','blue','green','black'};
        colororder(newcolors)
        fig1 = semilogy(ps,detP_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        fig1 = semilogy(ps,detP_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        % fig1 = semilogy(ps,detP_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
        % hold on
        fig1 = semilogy(ps,detP_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
        hold off
        legend({}, 'FontSize', 14)
        ax = gca; ax.FontSize = 14;
        xlabel('Number of sensors', 'FontSize', 14)
        ylabel('det P', 'FontSize', 14)
        filename = [ 'loop', num2str(w), '_detP.png' ];
        saveas(fig1, filename);
        
        colororder(newcolors)
        fig2 = semilogy(ps,detWo_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        fig2 = semilogy(ps,detWo_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        % fig1 = semilogy(ps,detWo_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
        % hold on
        fig2 = semilogy(ps,detWo_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
        hold off
        legend({}, 'FontSize', 14, 'Location', 'southeast')
        ax = gca; ax.FontSize = 14;
        xlabel('Number of sensors', 'FontSize', 14)
        ylabel('det W_o', 'FontSize', 14)
        filename = [ 'loop', num2str(w), '_detWo.png' ];
        saveas(fig2, filename);
        
        colororder(newcolors)
        fig3 = semilogy(ps,detCC_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        fig3 = semilogy(ps,detCC_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        % fig1 = semilogy(ps,detCC_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
        % hold on
        fig3 = semilogy(ps,detCC_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
        hold off
        legend({}, 'FontSize', 14, 'Location', 'southeast')
        ax = gca; ax.FontSize = 14;
        xlabel('Number of sensors', 'FontSize', 14)
        ylabel('det CC^T or det C^TC', 'FontSize', 14)
        filename = [ 'loop', num2str(w), '_detCC.png' ];
        saveas(fig3, filename);
        
        colororder(newcolors)
        fig4 = semilogy(ps,1./detP_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        fig4 = semilogy(ps,1./detP_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
        hold on
        % fig1 = semilogy(ps,1./detP_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
        % hold on
        fig4 = semilogy(ps,1./detP_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
        hold off
        legend({}, 'FontSize', 14, 'Location', 'southeast')
        ax = gca; ax.FontSize = 14;
        xlabel('Number of sensors', 'FontSize', 14)
        ylabel('det P^{-1}', 'FontSize', 14)
        filename = [ 'loop', num2str(w), '_detP-1.png' ];
        saveas(fig4, filename);
        cd(srcdir)
    end
    
    %% Data organization ===============================
    % Averaging
    [ time_KF, detWo_KF, detWo_KF_2, detWo_KF_3, detP_KF, detCC_KF ]...
        = F_data_ave1( CNT, num_ave, time_KF, detWo_KF, detWo_KF_2, detWo_KF_3, detP_KF, detCC_KF );
    [ time_std_KF, detWo_std_KF, detWo_std_KF_2, detWo_std_KF_3, detP_std_KF, detCC_std_KF ]...
        = F_data_std1( CNT, num_ave, time_KF, detWo_KF, detWo_KF_2, detWo_KF_3, detP_KF, detCC_KF );
    [ time_Gram, detWo_Gram, detWo_Gram_2, detWo_Gram_3, detP_Gram, detCC_Gram ]...
        = F_data_ave1( CNT, num_ave, time_Gram, detWo_Gram, detWo_Gram_2, detWo_Gram_3, detP_Gram, detCC_Gram );
    [ time_std_Gram, detWo_std_Gram, detWo_std_Gram_2, detWo_std_Gram_3, detP_std_Gram, detCC_std_Gram ]...
        = F_data_std1( CNT, num_ave, time_Gram, detWo_Gram, detWo_Gram_2, detWo_Gram_3, detP_Gram, detCC_Gram );
    [ time_S2SwR, detWo_S2SwR, detWo_S2SwR_2, detWo_S2SwR_3, detP_S2SwR, detCC_S2SwR ]...
        = F_data_ave1( CNT, num_ave, time_S2SwR, detWo_S2SwR, detWo_S2SwR_2, detWo_S2SwR_3, detP_S2SwR, detCC_S2SwR );
    [ time_std_S2SwR, detWo_std_S2SwR, detWo_std_S2SwR_2, detWo_std_S2SwR_3, detP_std_S2SwR, detCC_std_S2SwR ]...
        = F_data_std1( CNT, num_ave, time_S2SwR, detWo_S2SwR, detWo_S2SwR_2, detWo_S2SwR_3, detP_S2SwR, detCC_S2SwR );
    [ time_S2S, detWo_S2S, detWo_S2S_2, detWo_S2S_3, detP_S2S, detCC_S2S ]...
        = F_data_ave1( CNT, num_ave, time_S2S, detWo_S2S, detWo_S2S_2, detWo_S2S_3, detP_S2S, detCC_S2S );
    [ time_std_S2S, detWo_std_S2S, detWo_std_S2S_2, detWo_std_S2S_3, detP_std_S2S, detCC_std_S2S ]...
        = F_data_std1( CNT, num_ave, time_S2S, detWo_S2S, detWo_S2S_2, detWo_S2S_3, detP_S2S, detCC_S2S );
    
    % Arrange
    [time_all]  = F_data_arrange1( ps, CNT, time_KF, time_Gram, time_S2SwR, time_S2S );
    [detWo_all] = F_data_arrange1( ps, CNT, detWo_KF, detWo_Gram, detWo_S2SwR, detWo_S2S );
    [detWo2_all]= F_data_arrange1( ps, CNT, detWo_KF_2, detWo_Gram_2, detWo_S2SwR_2, detWo_S2S_2 );
    [detWo3_all]= F_data_arrange1( ps, CNT, detWo_KF_3, detWo_Gram_3, detWo_S2SwR_3, detWo_S2S_3 );
    [detP_all]  = F_data_arrange1( ps, CNT, detP_KF, detP_Gram, detP_S2SwR, detP_S2S );
    [detCC_all] = F_data_arrange1( ps, CNT, detCC_KF, detCC_Gram, detCC_S2SwR, detCC_S2S );
    [time_std]  = F_data_arrange1( ps, CNT, time_std_KF, time_std_Gram, time_std_S2SwR, time_std_S2S );
    [detWo_std] = F_data_arrange1( ps, CNT, detWo_std_KF, detWo_std_Gram, detWo_std_S2SwR, detWo_std_S2S );
    [detWo2_std]= F_data_arrange1( ps, CNT, detWo_std_KF_2, detWo_std_Gram_2, detWo_std_S2SwR_2, detWo_std_S2S_2 );
    [detWo3_std]= F_data_arrange1( ps, CNT, detWo_std_KF_3, detWo_std_Gram_3, detWo_std_S2SwR_3, detWo_std_S2S_3 );
    [detP_std]  = F_data_arrange1( ps, CNT, detP_std_KF, detP_std_Gram, detP_std_S2SwR, detP_std_S2S );
    [detCC_std] = F_data_arrange1( ps, CNT, detCC_std_KF, detCC_std_Gram, detCC_std_S2SwR, detCC_std_S2S );
    
    
    %% Save =============================================================
    cd(workdir)
    % Average
    save('time.mat','time_all');
    save('detWo.mat','detWo_all');
    save('detWo2.mat','detWo2_all');
    save('detWo3.mat','detWo3_all');
    save('detP.mat','detP_all');
    save('detCC.mat','detCC_all');
    % std
    save('time_std.mat','time_std');
    save('detWo_std.mat','detWo_std');
    save('detWo2_std.mat','detWo2_std');
    save('detWo3_std.mat','detWo3_std');
    save('detP_std.mat','detP_std');
    save('detCC_std.mat','detCC_std');
    % each detP
    save('detP_KF.mat','detP_KF');
    save('detP_Gram.mat','detP_Gram');
    save('detP_S2SwR.mat','detP_S2SwR');
    save('detP_S2S.mat','detP_S2S');
    % each detWo
    save('detWo_KF.mat','detWo_KF');
    save('detWo_Gram.mat','detWo_Gram');
    save('detWo_S2SwR.mat','detWo_S2SwR');
    save('detWo_S2S.mat','detWo_S2S');
    % each detWo2
    save('detWo_KF_2.mat','detWo_KF_2');
    save('detWo_Gram_2.mat','detWo_Gram_2');
    save('detWo_S2SwR_2.mat','detWo_S2SwR_2');
    save('detWo_S2S_2.mat','detWo_S2S_2');
    % each detWo3
    save('detWo_KF_3.mat','detWo_KF_3');
    save('detWo_Gram_3.mat','detWo_Gram_3');
    save('detWo_S2SwR_3.mat','detWo_S2SwR_3');
    save('detWo_S2S_3.mat','detWo_S2S_3');
    % each detCC
    save('detCC_KF.mat','detCC_KF');
    save('detCC_Gram.mat','detCC_Gram');
    save('detCC_S2SwR.mat','detCC_S2SwR');
    save('detCC_S2S.mat','detCC_S2S');
    % each time
    save('time_KF','time_KF');
    save('time_Gram','time_Gram');
    save('time_S2SwR','time_S2SwR');
    save('time_S2S','time_S2S');
        
    diary off
    tEnd=toc(tStart)
    save('runtime.mat','tEnd')
    cd(srcdir)
    disp('All results are saved!');
    

%% Practical problem using PIV dataset ==========================
elseif num_problem == 3
    % Loading PIV data
    text='Readinng/Arranging a PIV dataset'; disp(text);
    load([datadir,'/',num2str(test_date),'_run',num2str(run_num),'_0_41A_80us_10ms_',num2str(AoA),'deg.mat']);
    % (this file includes: AUN, AveU, AveW, AWN, IMask_U, Nend, Nstart, Uall, Wall)

    % Generation of data matrix except mask area (copy Kanda-code）
    [a,b,nT] = size(Uall);
    g = 0;
    for e = 1:a
        for f = 1:b
            if IMask_U(e,f) == 1
                g = g + 1;
            end
        end
    end
    Xorg = zeros(2*g,nT);
    xv = zeros(g,nT);
    xv_temp = zeros(a*b,1);
    yv = zeros(g,nT);
    yv_temp = zeros(a*b,1);
    xave = zeros(g,nT);
    xave_temp = zeros(a*b,1);
    yave = zeros(g,nT);
    yave_temp = zeros(a*b,1);
    IMask_vec = reshape(IMask_U(:,:),[],1);
    for i = 1:nT
        xv_temp = reshape(Uall(:,:,i),[],1);
        yv_temp = reshape(Wall(:,:,i),[],1);
        xave_temp = reshape(AveU,[],1);
        yave_temp = reshape(AveW,[],1);
        %ProgressBar(i,nT)
        k = 1;
        for j = 1: a*b
            if IMask_vec(j) == 1
                xv(k,i) = xv_temp(j);
                yv(k,i) = yv_temp(j);
                xave(k,i) = xave_temp(j);
                yave(k,i) = yave_temp(j);
                k = k + 1;
            end
        end

    end
    Xorg(1:g, :) = xv - xave;
    Xorg((g+1):(2*g), :) = yv - yave;

    % POD (copy and modify Kanda-code)
    [U,S,V] = svd(Xorg);
    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Vr = V(:,1:r);
    nall = size(Xorg,1);%nall is defined as product of n and s(2021/11/15)
    n    = nall/s;
    m    = size(Xorg,2);
    save([workdir2,'/X.mat'],'Xorg')%, 'xave', 'IMask_vec')
    save([workdir2,'/U.mat'],'U')
    save([workdir2,'/S.mat'],'S')
    save([workdir2,'/V.mat'],'V')
    save([workdir2,'/xave.mat'],'xave')
    save([workdir2,'/yave.mat'],'yave')
    save([workdir2,'/IMask_vec.mat'],'IMask_vec')

    % Generation of Aorg and Corg
    Z    = Sr*Vr';
    Aorg = Z(:,2:end) * pinv(Z(:,1:end-1)); %Kanda-code F
    Corg = Ur;        %case of all s
    Borg = zeros(r,1);
    
    % - Case1: Load Q,R calculated previously using P0_estimateQR -----
    text='Loading QR matrices previously calculated'; disp(text);
    load([workdir3,'/Q.mat']);
    load([workdir3,'/R.mat']);
    
    % - Case2: Calculate Q,R using function(Note that this takes long time)
%     text='Calculating QR matrices using function'; disp(text);
%     [Q,R] = F_estimate_QR(Xorg,Ur,Sr,Vr,Aorg,workdir3);
%     [Q,R] = F_estimate_QRdense(Xorg,Ur,Sr,Vr,Aorg,workdir3);

    %% Sensor selection w/o average loop
    p = pmax;
    
    % Dynamics w/ noise -------------------------------
    % - diag(QR), Vector
    % [time_KF, H_KF, sensors_KF, ranki] = F_sensor_KF_QRvec(p,s,Aorg,Borg,Corg,Q,R);
    % - dense(QR), Vector
    [time_KF, H_KF, sensors_KF, ranki] = F_sensor_KF_denseQRvec(p,s,Aorg,Borg,Corg,Q,R);
    
    % S2S w/o noise -----------------------------------
    % Maximization of determinant of Fisher information matrix
    % - Vector
    [time_S2S, H_S2S, sensors_S2S, det_test] = F_sensor_DGvec(p,s,Corg);
    
    %% Evaluation of optimality indices
    CNT = 0;
    text = [ 'Evaluation started --->' ]; disp(text);
    for p=ps
        CNT=CNT+1;
        detP_KF   (CNT,1) = F_calc_detP_QRvec(p,s,H_KF(1:p,1:n),sensors_KF(1:p),Aorg,Borg,Corg,Q,R);
        detCC_KF  (CNT,1) = F_calc_det_vec   (p,s,H_KF(1:p,1:n),Corg);
        % detWo code is still not developed...
        
        detP_S2S   (CNT,1) = F_calc_detP_QRvec(p,s,H_S2S(1:p,1:n),sensors_S2S(1:p),Aorg,Borg,Corg,Q,R);
        detCC_S2S  (CNT,1) = F_calc_det_vec   (p,s,H_S2S(1:p,1:n),Corg);
        
    end
    text = [ '---> Evaluation finished.' ]; disp(text);
    
    %% Make R matrix
    [R_KF]  = F_make_R(R,s,sensors_KF);
    [R_S2S] = F_make_R(R,s,sensors_S2S);
    
    %% Data organization
    % Arrange
    [detP_all]  = F_data_arrange2( ps, CNT, detP_KF, detP_S2S );%detP_Gram, detP_S2SwR,
    [detCC_all] = F_data_arrange2( ps, CNT, detCC_KF, detCC_S2S );%detCC_Gram, detCC_S2SwR,
    
    %% Save
    cd(workdir3)
    cond = ['p',num2str(pmax)];
    save(['H_KF_',cond,'.mat'],'H_KF')
    save(['sensors_KF_',cond,'.mat'],'sensors_KF')
    save(['R_KF_',cond,'.mat'],'R_KF')
    save(['H_S2S_',cond,'.mat'],'H_S2S')
    save(['sensors_S2S_',cond,'.mat'],'sensors_S2S')
    save(['R_S2S_',cond,'.mat'],'R_S2S')
    % all
    save(['detP_',cond,'.mat'],'detP_all');
    save(['detCC_',cond,'.mat'],'detCC_all');
    % each detP
    save(['detP_KF_',cond,'.mat'],'detP_KF');
    save(['detP_S2S_',cond,'.mat'],'detP_S2S');
    % each detCC
    save(['detCC_KF_',cond,'.mat'],'detCC_KF');
    save(['detCC_S2S_',cond,'.mat'],'detCC_S2S');
    cd(srcdir)
    disp('All results are saved!')
    
end    


%% /////////////////////////////////////////////////////////////
% Main program end


