% function to select sensor based on KF
% F_sensor_KF.m : ------------- Programer : Kumi Nakai, Taku Nonomura 2020
%                               Last modified: 2021/11/17 K.Nakai(maintanence)
% Noise (system&observation): White noise is produced using amplitude given as input param
% Component number: 1

function [time, isensors]=F_sensor_WSDP_KF(Aorg,Corg,Q,sigma_s2,p)%,Dorg
n=size(Corg,1);
tic
[zhat,z,X,Y,Z,sensors]=F_sensor_WSDP_KF_CVX(Aorg,Corg,Q,sigma_s2,p);
isensors = sensors(1:p);
%[H]=F_calc_sensormatrix(p, n, sensors);
D=0;
time=toc;  

end

