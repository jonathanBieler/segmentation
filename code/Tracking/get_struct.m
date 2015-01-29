function mystruct = get_struct(type, nstruct)
% GET_STRUCT retrieve custom data structures.
%   This function is designed as a centralization for the different 
%   complex structures used throughout the program so that they can
%   be edited consistently.
%
%   MYSTRUCT = GET_STRUCT(TYPE, SIZE) returns a matrix of size SIZE 
%   of the custom data structure of type TYPE. SIZE can be multi-
%   dimensional.
%
%   MYSTRUCT = GET_STRUCT(TYPE) returns one structure (SIZE = 1).
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 09.12.2010

  % Set the default size
  if (nargin == 1)
    nstruct = 1;
  end

  % Switch between all the different types
  switch type

    % The general structure that contains all the information required for
    % ASSET to analyze the data. Usually this structure is called 'opts' throughout the code.
    case 'ASSET'
      mystruct = struct('analyzed_fields', {{'carth'}}, ... % Fields in mymovie used for the analyis
                        'application', {{''}}, ...          % List of the applications (other than segmentation) that will be performed
                        'auto_save', true, ...              % Automatically save the intermediate results
                        'binning', 1, ...                   % Pixel binning used during acquisition
                        'ccd_pixel_size', 1, ...         % X-Y size of the pixels in ?m (of the CCD camera, without magnification)
                        'compression', 'none', ...           % Compression used for the data files (prompted is empty)
                        'compute_probabilities', false, ... % Compute the posterior probability
                        'config_file', '', ...              % Name of the configuration file that will be loaded 
                        'crop_export', false, ...           % Crop the images when exporting the results
                        'crop_size', 2.2, ...               % Crop size is the axes_length*crop_size
                        'debug', false, ...                 % Debug mode ON
                        'do_ml', 'none', ...                % Machine learning (ML) is performed
                        'dp_method', 'double', ...          % Dynamic programming method used (see dynamic_programming.m)
                        'export_movie', false, ...          % Export the results of the analysis
                        'file_regexpr', '(.+[-_])?(.*?[-_]?)(\.[\w\.]+)', ... % The regular expression representing the files
                        'filters', get_struct('channel_filter'), ...
                        'follow_periphery', true, ...       % Follow only the periphery of the trackings (no invaginations)
                        'force_circularity', true, ...      % Enforces that the last row of the DP conincides with the first one
                        'magnification', 1, ...            % Magnification of the objective of the microscope
                        'max_export', 1, ...                % Number of exported frames (>1 is the absolute number of frames, <=1 is the fraction)
                        'max_frames', 1, ...                % Number of analyzed frames (value as for max_export)
                        'measure_performances', false, ...  % Measures the error of the segmentation with respect to the manual trackings
                        'merge_input_files', true, ...      % Uses the Bio-formats file merge when importing movie files
                        'ml_type', 'cortex', ...            % Field on which ML is performed (use cell array for more than one)
                        'nbins', 36, ...                    % Number of bins used to measure the performance
                        'normalize', true, ...              % Normalize the results of the analysis onto the reference embryo
                        'overwrite', true, ...              % Overwrite the previous data by saving in the same MAT-file
                        'parse_export', 'normal', ...       % How export is performed (normal or random)
                        'parse_frames', 'normal', ...       % Order of the frames for the segmentation (normal or random)
                        'pixel_size', 0, ...                % X-Y size of the pixels in ???m (computed as ccd_pixel_size / magnification)                        
                        'recompute', false, ...             % Recompute previously computed features (mainly segmentaiton and trackings)                       
                        'spot_tracking', get_struct('spot_tracking'), ... % Parameters of the spot tracking algorithm                        
                        'trackings', '', ...                % List of tracking files
                        'uuid', 0 , ...                     % Universal Unique IDentifier (used in ML to identify processes) 
                        'verbosity', 2, ...                 % Verbosity level (0 null, 1 text only, 2 gui, 3 full with plots)
                        'warp_type', 'radial');             % Warp type used to normalize the embryo (see carth2normalized.m)

      % Compute the pixel size based on the default values. This needs to be re-done
      % in case one value is changed.
      mystruct = set_pixel_size(mystruct);

    % Structure containing the different parameters required for tracking spots
    case 'spot_tracking'
      mystruct = struct('fusion_thresh', 2, ...         % Minimal distance in um to another spot (estimation) before fusion
                        'frame_displacement', 60, ...   % Maximal displacement of a spot (in random units) between two frames 
                        'split_cost',1,...              % Cost of splits: low cost -> more divisions
                        'frame_window', 2, ...          % Considered number of frames for the gap closing algorithm (see track_spots.m)
                        'gap_function', @relative_distance, ... % Function used to measure the gap-closing weight
                        'joining_function', @merging_distance, ... % Same but for the joinging weight
                        'splitting_function', @splitting_distance, ... % For the splitting weight
                        'linking_function', @mutual_distance, ... % And for the frame-to-frame linking 
                        'max_size', 30, ...              % Maximal size (in um) of the spots
                        'noise_thresh', 3);             % Threshold used to remove the nosie (see imatrou.m)

   
    % If the required type of structure has not been implemented, return an empty one
    otherwise
      mystruct = struct();
  end

  % Repeat the structure to fit the size (nstruct can be multi-dimensional)
  mystruct = repmat(mystruct,nstruct);

  return;
end
