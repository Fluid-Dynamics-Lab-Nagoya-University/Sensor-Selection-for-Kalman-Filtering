% function to select sensor based on KF
% F_sensor_KF.m : ------------- Programer : Kumi Nakai, Taku Nonomura 2020
%                               Last modified: 2021/11/17 K.Nakai(maintanence)
% Noise (system&observation): White noise is produced using amplitude given as input param
% Component number: 1

function [time, isensors]=F_sensor_GSDP_KF(Aorg,Corg,Q,sigma_s2,p)%,Dorg
n=size(Corg,1);
%lambda=0.5;  %30-200
%lambda=0.8; %300-3000

%lambda=0.1; % for1-30 sigma+3
%lambda=0.05; % for1-30 sigma+3
%lambda=0.5; % for1-30 sigma+1
%lambda=0.5; % for1-30 sigma+0
%lambda=0.5; % for1-30 sigma-1
%lambda=0.01;  % for1-30 sigma-3
%lambda=0.0001;  % for1-30 sigma-5

%sigma-1 p=20
%lambda = 0.15; %n=50
%lambda = 0.3; %n=100
%lambda = 0.4; %n=300
%lambda = 0.5; %n=500
%lambda = 0.6; %n=1000
%lambda = 0.5; %n=3000
lambda = 0.5; %n=5000

ppre=0;

while ppre < p
    isensors=[];
    tic
    [zhat,z,X,Y,Z,sensors]=F_sensor_GSDP_KF_CVX(Aorg,Corg,Q,sigma_s2,p,lambda);
    isensors = sensors(1:p);
    %[H]=F_calc_sensormatrix(p, n, sensors);
    D=0;
    time=toc;
    ppre=size(find(z>0.0001),1)
    if ppre < p
        lambda=lambda*0.8;
        disp('Retry with smaller lambda!')
    end
    disp(['ppre=',num2str(ppre)])
    disp(['lambda=',num2str(lambda)])
    pause(0.5)
end
end