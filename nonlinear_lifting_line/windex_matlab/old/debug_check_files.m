function check_files()
    files = {'wxwing_d-10_re3e5.txt', 'wxwing_d10_re3e5.txt', 'wxwing_re3e5.txt'};
    dir_path = '../windex_xfoil/';
    
    for i = 1:length(files)
        fname = fullfile(dir_path, files{i});
        if exist(fname, 'file')
            try
                raw = importdata(fname, ' ', 12);
                if isstruct(raw)
                    d = raw.data;
                else
                    d = readmatrix(fname, 'NumHeaderLines', 12);
                end
                alpha = d(:, 1);
                fprintf('%s: Alpha range [%.2f, %.2f], %d points\n', files{i}, min(alpha), max(alpha), length(alpha));
            catch
                fprintf('%s: Error reading\n', files{i});
            end
        else
            fprintf('%s: Not found\n', files{i});
        end
    end
end
