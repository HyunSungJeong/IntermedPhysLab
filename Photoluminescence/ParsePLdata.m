function PLdata = ParsePLdata(varargin)
    % <Description>
    % Loads photoluminescence data from .asc file
    %
    % <Option>
    % 'rhodamine' : If used, loads PL data for rhodamine at room temperature
    %                   (Default : used)
    %
    % 'ruby', ... : [numeric] Temperature at which the PL data for ruby to be loaded
    %                       If used, PL data for ruby is loaded
    %                   (Default: not used)
    %
    % 'rubyRtemp' : If used, PL data for ruby at room temperature is loaded
    %                   (Default: not used)
    %
    % 'photonE' : If used, the first column of the data is the single
    %               photon in meV unit, instead of wavelength in nm unit
    %                   (Default: not used)
    %
    % <Output>
    % PLdata : [numeric] Nx2 numeric array. 
    %                   The first column is wavelength in nm unit or energy of a single photon in meV unit,
    %                   and the second column is the CCD signal at each wavelength/energy

    %% Parse options

    % default options
    dataType = 'rhodamine';
    photonE = false;

    while ~isempty(varargin)
        switch varargin{1}

            case 'rhodamine'
                dataType = 'rhodamine';
                varargin(1) = [];

            case 'ruby'
                if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                    error('ERR: Temperature must be single positive number');
                elseif varargin{2} < 0
                    error('ERR: Absolute temperature cannot be negative');
                else
                    dataType = 'ruby';
                    T = varargin{2};
                    varargin(1:2) = [];
                end

            case 'rubyRtemp'
                dataType = 'rubyRtemp';
                varargin(1) = [];

            case 'photonE'
                photonE = true;
                varargin(1) = [];

            otherwise
                if ischar(varargin{1})
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                else
                    error('ERR: Unknown input');
                end
        end
    end

    %% Read file

    if isequal(dataType, 'rhodamine')
        filename = [fileparts(mfilename('fullpath')), filesep, 'PLdata', filesep, 'rhodamine.asc'];
    elseif isequal(dataType, 'ruby')
        filename = [fileparts(mfilename('fullpath')), filesep, 'PLdata', filesep, 'ruby_', sprintf('%d', T), 'K.asc'];
    else
        filename = [fileparts(mfilename('fullpath')), filesep, 'PLdata', filesep, 'ruby_roomtemp.asc'];
    end
    
    % Open the file
    fid = fopen(filename, 'r');
    
    % Check if file opened successfully
    if fid == -1
        error('Cannot open the file: %s', filename);
    end
    
    % Skip lines until we find the start of numeric data
    while ~feof(fid)
        pos = ftell(fid); % Remember current file position
        line = fgetl(fid);
        
        % Try parsing the line as two numeric values
        numbers = sscanf(line, '%f %f');
        
        % If successful (i.e., we get two values), rewind and break
        if length(numbers) == 2
            fseek(fid, pos, 'bof'); % rewind to the beginning of the line
            break;
        end
    end
    
    % Now read the remaining data (assumes two columns of numeric data)
    PLdata = fscanf(fid, '%f %f', [2 Inf]);
    
    % Close the file
    fclose(fid);
    
    % Transpose data to get two columns
    PLdata = PLdata';

    if photonE
        h = 6.626 * 1e-34;
        c = 3 * 1e8;
        e_charge = 1.602 * 1e-19;
        PLdata(:,1) = 1e12 * h*c./PLdata(:,1) .* (1/e_charge);
    end

end