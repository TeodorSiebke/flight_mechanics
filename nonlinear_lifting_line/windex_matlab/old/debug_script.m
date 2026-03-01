addpath('../matlab');
delete('debug_output.txt');
diary('debug_output.txt');

try
    acg = make_windex();
    uoo = 30;
    fs = flight_state(acg, uoo, 0, 0, 0);
    
    fprintf('Fields of fs:\n');
    fnames = fieldnames(fs);
    for i=1:length(fnames)
        fprintf('  %s\n', fnames{i});
    end
    
    alpha_deg = (-2:1:12)';
    alpha_rad = alpha_deg * pi/180;
    
    fprintf('\nStarting plainpolar...\n');
    [CL, CD, CM] = plainpolar(acg, fs, alpha_rad);
    
    fprintf('\nResults alpha, CL, CM:\n');
    for i=1:length(alpha_deg)
        fprintf('data_point: %.1f, %.4f, %.4f\n', alpha_deg(i), CL(i), CM(i));
    end
    
catch ME
    fprintf('\nERROR: %s\n', ME.message);
    for i=1:length(ME.stack)
        fprintf('  in %s at line %d\n', ME.stack(i).name, ME.stack(i).line);
    end
end

diary off;
exit;
