function temp = get_tctl_temperature(folder_path)
% GET_TCTL_TEMPERATURE Retrieves the GPU's Tctl temperature using the 'sensors' command.
%
%   TEMP = GET_TCTL_TEMPERATURE() runs the 'sensors' command, parses its output
%   for the line containing "Tctl:" and extracts the temperature in °C.
%
%   TEMP = GET_TCTL_TEMPERATURE(FOLDER_PATH) does the same, but if FOLDER_PATH is provided
%   and not empty, it logs the current temperature along with the current date and time
%   to a log file named "temp_log" in the specified folder. If the file exists, the new
%   data is appended.
%
%   If the temperature is found, the function returns it as a double.
%   Otherwise, it returns NaN.

    % Check if folder_path is provided; if not, set to empty.
    if nargin < 1
        folder_path = '';
    end

    try
        % Execute the 'sensors' command
        [status, output] = system('sensors');
        if status ~= 0
            error('Error running sensors command');
        end

        % Define the regular expression to capture the temperature.
        % Pattern explanation:
        %   'Tctl:'   - looks for the literal string "Tctl:"
        %   '\s+'     - one or more whitespace characters
        %   '\+'      - a literal plus sign
        %   '([\d.]+)'- one or more digits or a decimal point (captured)
        %   '°C'      - literal degree symbol followed by C
        expr = 'Tctl:\s+\+([\d.]+)°C';
        tokens = regexp(output, expr, 'tokens');

        if ~isempty(tokens)
            % Extract the first match
            tempStr = tokens{1}{1};
            temp = str2double(tempStr);
        else
            % No match found
            temp = NaN;
        end

    catch ME
        warning('getTctlTemperature:tempRetrieveError', ...
            'Error retrieving temperature: %s', ME.message);
        temp = NaN;
    end

    % If a folder path is provided, log the current temperature and time.
    if ~isempty(folder_path)
        % Create the full path for the log file
        log_filename = fullfile(folder_path, 'temp_log');
        
        % Get the current timestamp in a readable format
        tstamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        
        % Prepare the log entry (temperature formatted to two decimal places)
        log_entry = sprintf('%s: %.2f°C\n', tstamp, temp);
        
        % Open the file in append mode ('a') and write the log entry
        fid = fopen(log_filename, 'a');
        if fid == -1
            warning('getTctlTemperature:logFileError', 'Could not open log file: %s', log_filename);
        else
            fprintf(fid, '%s', log_entry);
            fclose(fid);
        end
    end

end
