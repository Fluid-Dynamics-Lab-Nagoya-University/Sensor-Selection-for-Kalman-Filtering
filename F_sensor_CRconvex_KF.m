function [time, isensors, data]=F_sensor_CRconvex_KF(A, U, Q, sigma2, p, maxiteration,l)
n=size(U,1);
tic
[zhat, ~, ztilde, ~, sensors, data ] = F_sens_sel_KF_CRapproxnt(A, U, Q, sigma2, p, maxiteration,l);
%maximum selection
isensors = sensors(1:p);
%[H]=F_calc_sensormatrix(p, n, sensors);
D=0;
time=toc;