function [time, sensors]=F_sensor_random(n,p)

    tic;
    sensors = randperm(n,p);
    time=toc;
    %[H]=F_calc_sensormatrix(p, n, sensors);
    
end