addpath('../matlab');
acg = make_windex();
acg.pref = [0.25, 0.0, 0.0];
uoo = 90;
fs = flight_state(acg, uoo, 0, 0, 0);
fs.delta_flap = 0;
fs.delta_aileron = 0; fs.delta_rudder = 0;
fs.delta_elevator = 0;

fprintf('Testing 90 m/s at Alpha=0...\n');
try
    gamma = pointsolve(acg, fs);
    fs.gamma = gamma; 
    [CL, CD, CY, Cm] = coefficients(acg, fs);
    fprintf('Success! CL=%.4f, CD=%.4f, Cm=%.4f\n', CL, CD, Cm);
catch ME
    fprintf('Failed: %s\n', ME.message);
    for k=1:length(ME.stack)
        fprintf('  in %s at %d\n', ME.stack(k).name, ME.stack(k).line);
    end
end
