function load_data()
    script_dir = fileparts(mfilename('fullpath'));
    
    files = {'wxwing_d-10_re3e5.txt', 'wxwing_d10_re3e5.txt', 'wxwing_re3e5.txt'};
    dir_path = fullfile(script_dir, '../windex_xfoil/');
    
    out_file = fullfile(script_dir, 'output_log.txt');
    fid = fopen(out_file, 'w');
    fprintf(fid, 'Starting load_data...\n');
    
    for i = 1:length(files)
        fname = fullfile(dir_path, files{i});
        if exist(fname, 'file')
            try
                % Import data
                d = importdata(fname, ' ', 12);
                if isstruct(d)
                    dat = d.data;
                else
                    dat = d;
                end
                
                % Check if dat is empty
                if isempty(dat)
                     fprintf(fid, '%s: Empty data\n', files{i});
                else
                     fprintf(fid, '%s: Range [%.2f, %.2f], %d points\n', ...
                        files{i}, min(dat(:,1)), max(dat(:,1)), length(dat(:,1)));
                end
            catch ME
                fprintf(fid, '%s: Error reading: %s\n', files{i}, ME.message);
            end
        else
            fprintf(fid, '%s: Not found in %s\n', files{i}, dir_path);
        end
    end
    fprintf(fid, 'Finished load_data.\n');
    fclose(fid);
end
