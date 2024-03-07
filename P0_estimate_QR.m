% P0_estimate_QR.mat ----------- Programer : Kumi Nakai 2021
%                               Last modified: 2021/11/18 K.Nakai
% Load training data and calculate Q,R matrices 
% using all observation points arranged in the order of initial condition.
% You can choose to run this code before running P1_dynamics_PIV or to use F_estimate_QR
% during running P1_dynamics_PIV.
% F_estimate_QR.m: P0_estimate_QR is functionalized to adapt to P1_dynamics_PIV.

clearvars
close all;

%% Input parameters
s =  2;
r = 10;
test_date = 20210518; run_num = 148; AoA = 18;

%% Preparation of output directories ================================
srcdir  = pwd;
datadir = '../data';
workdir = ['../work_sensor/',num2str(AoA),'deg_s',num2str(s),'/r',num2str(r)];
if not(exist(workdir,'dir'))
    mkdir(workdir); disp(['New directory "',workdir,'" was made']);
end

%% Load training data (Not necessary in Function(F_estimateQR.m))
% Kanda-code ===
% filename = ([num2str(test_date),'_run',num2str(run_num),'_0_41A_80us_10ms_',num2str(AoA),'deg_sensorDG']);
% load([filename, 'U.mat']);
% load([filename, 'S.mat']);
% load([filename, 'V.mat']);
% %load(['20200722_run19_0_41A_10ms_14deg_sensor_',num2str(p),'_DG_H.mat']);
% if r == 10
%     load([filename, '_sensor_',num2str(p),'H.mat']);
% else
%     load([filename,'_r_',num2str(r), '_sensor_',num2str(p),'H.mat']);
% end
% dirname = ('Analyze_data\Kalman_data_fft2');
% mkdir(dirname);
% dirname = (['Analyze_data\Kalman_data_fft2\',filename]);
% mkdir(dirname);
% ===

cd(datadir)
filename = [num2str(test_date),'_run',num2str(run_num),'_0_41A_80us_10ms_',num2str(AoA),'deg'];
load([filename,'.mat']);
cd(srcdir)%('../src')

% (Copy Kanda-code)
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

% Added by Nakai 2021/11/15
nall = size(Xorg,1); %nall is defined as product of n and s(2021/11/15)
n    = nall/s;
m    = size(Xorg,2);

% POD (Copy and modify Kanda-code)
[U,S,V] = svd(Xorg);
Ue = U(:,1:r);
Se = S(1:r,1:r);
Ve = V(:,1:r);
Z = Se*Ve';
F = Z(:,2:end) * pinv(Z(:,1:end-1));
C = Ue;%H * Ue; !!!

%% P,Q,R calculation (Copy and modify Kanda-code)
P = eye(r,r);
Q = eye(r,r);
R = eye(nall,nall);

for mm = 1:m
    y = Xorg(:,mm);%H * Xorg(:,mm) !!!
    z = Z(:,mm);
    P = F*P*F'+Q;
    K = P*C'* inv(C*P*C' + R);
    P = P - K*C*P;
   
    [R]=F_estimate_R_sub(C, z, y, nall);
    Restimate(1:nall,nall*(mm-1)+1:nall*(mm-1)+nall)=R;
    if mm < m
        [Q]=F_estimate_Q_sub(F, Z, r, mm, m);
        Qestimate(1:r,r*(mm-1)+1:r*(mm-1)+r)=Q;
    end
    Q = eye(r,r);
    R = eye(nall,nall);
end

[Q,R] = F_estimate_QR_sub(Qestimate, Restimate, nall, r, m);

save([workdir,'/Q.mat'],'Q')
save([workdir,'/R.mat'],'R')

