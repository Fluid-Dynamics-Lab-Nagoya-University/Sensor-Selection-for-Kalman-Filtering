function [time, isensors]=F_sensor_convex_KF(A, U, Q, sigma2, p, maxiteration)
n=size(U,1);
tic
[zhat, ~, ztilde, ~, sensors, data ] = F_sens_sel_KF_approxnt(A, U, Q, sigma2, p, maxiteration);
%maximum selection
isensors = sensors(1:p);
%[H]=F_calc_sensormatrix(p, n, sensors);
D=0;
time=toc;