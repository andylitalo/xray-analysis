% Use this script if you would like to process several scans and average
% all of the frames within each scan of interest.

% delete if MATLAB starts in your current directory--otherwise, put the
% file path to your current director in the string after "cd"
% cd <current directory here>

% clear MATLAB
clear all
close all
clc

% FILE PARAMETERS
% file name until scan number
file_hdr = 'eimear_hs104_'; 
% image extension
ext = '.tif';
% number of the scans that you want to analyze
scan_arr = [236, 237];
% number of frames per scan
n_frames = 3;
% average signal frames from a given scan?
avg_frames = true;
% number of the scan used for background, typically the glass mold scan
bkg_scan = 235;
% number of frames in the background scan
n_bkgd = 5;
% max intensity of images to scale 2D bkgd-subt images [a.u.]
% change this to adjust the color scale on the plot
max_intensity = 1000;

% FIGURE SETTINGS
% detector type (SAXS, MAXS, or WAXS)--just affects naming of files & plots
detector_type = 'SAXS';
% save the plot of the 1D integrated intensity?
save_I_1D = true;
% open a new window for each folder?
open_new_windows = true;
% two parameter subtraction?
two_para_subt = false;

% (GENERALLY) FIXED PARAMETERS
% analysis settings
% min and max radii from beam stop to analyze [pixels]
if strcmp(detector_type, 'SAXS')
    r_min = 10;
    r_max = 200; 
    % min and max azimuthal angles to analyze [degrees]
    phi_min = -180;
    phi_max = 180;
elseif strcmp(detector_type, 'WAXS')
    r_min = 0;
    r_max = 1000;
    % min and max azimuthal angles to analyze [degrees]
    phi_min = 0;
    phi_max = 180;
else
    fprintf('Define parameters for detector type %s.\n', detector_type)
end
% Intervals of the wave number q for 2 parameter subtraction [1/A] 
% these should be regions where the background and signal agree
q_low_min = 0.70; 
q_low_max = 0.75; 
q_high_min = 1.80; 
q_high_max = 1.90; 
% original beam center coordinates [pixels]--you can find this in the spec
% file
row_center = 975; % y center
col_center = 960; % x center
% Beam parameters
lamda = 0.7293*10^-10; % wavelength [m]
pixsize = 88.6 *10^-6; % pixel size [m]  
SaDet = 200.51 * 10^-3 ; % sample to detector distance [m]
exposure_time = 0.5; % exposure/acquisition time [s]
parameters = [lamda; pixsize; SaDet];
% font parameters
set(0,'defaultAxesFontSize',20);
set(0,'defaultTextFontSize',20);
set(0,'defaultTextFontName','Times');
% font size of legend
l_fs = 12;

% directories
dir_1d = '..\\1d_integrated\\';
dir_2d_bkgd_subt = '..\\2d_bkgd_subt\\';
dir_2d = '..\\2d_raw\\';
dir_data = '..\\datasheets\\';

%% set parameters for integration
 
% find size of image
filename = strcat(dir_2d, file_hdr, num2str(scan_arr(1)), '_0000', ext); 
[im_row, im_col] = size(imread(filename));
% % Parameters for radius and phi matrices
nr = r_max - r_min; % number of points along the radial direction
nphi = round (2*pi*r_max); % number of points along the azimuthal direction

% create a rectangular matrix of coordinates for phi and the radius
[phi, radius] = meshgrid (linspace(phi_min, phi_max, nphi), linspace(r_min,r_max,nr));
% create a rectangular matrix of x, y coordinates of the same size 
[x,y] = meshgrid(1:max(im_col, col_center),1:max(im_row, row_center)); 

% convert polar to cartesian coordinates
% Note that r*sin(phi) and r*cos(phi) are subtracted instead of added
% because phi is measured from the left, contrary to convention from the
% right
yi = row_center - radius.*sind(phi);
xi = col_center - radius.*cosd(phi);

%% Average the scans of the background (typically the glass mold)
bkg_sum = 0; 
% THIS WILL FAIL IF YOUR BACKGROUNDS ARE NOT INDEXED FROM 0000
for k = 0:n_bkgd-1
    % read file and append rows
    filename = strcat(dir_2d, file_hdr, num2str(bkg_scan), sprintf('_%04d', k), ext); 
    rows_needed = row_center - im_row; 
    % load pixel values of image--THIS WILL FAIL IF YOU NEED TO ADD COLUMNS
    img = vertcat(double(imread(filename)), zeros(rows_needed, im_col)); 
    % add background scan
    bkg_sum = bkg_sum + img; 
end
% divide by number of background scans to get the average
bkg_avg = bkg_sum./n_bkgd; 

%% Subtract background from samples; save 2D and 1D data 

if avg_frames
    n_frames_2_save = 1;
else
    n_frames_2_save = n_frames;
end
% store data with background subtracted
IvsQ_Sub = zeros(nr,length(scan_arr)*n_frames_2_save);
IvsPhi_Sub = zeros(nphi,length(scan_arr)*n_frames_2_save);
% store data without subtracting background, plus store the background
IvsQ = zeros(nr, length(scan_arr)*n_frames_2_save+1);
IvsPhi = zeros(nphi, length(scan_arr)*n_frames_2_save+1);
% Initialize data structures
labels = {}; % labels for legend
alpha = zeros(length(scan_arr)*n_frames_2_save,1); % parameter 1 for 2-param subtraction
beta = zeros(length(scan_arr)*n_frames_2_save,1); % parameter 2 for 2-param subtraction
% calculate the q-matrix (for calibration from x-y to q-phi space in
% 2-parameter subtraction)
q_matrix = Q_Matrix_Calculator();
% initialize header for data files
hdr = zeros(1, length(scan_arr)*n_frames_2_save);
% Analyze samples
% (try parallel for-loop -> this uses all cores for non-sequential calcs)
ct = 1;
for k = 1:length(scan_arr)
    % Load sample number
    scan = scan_arr(k); 
    disp (scan) 
    
    if avg_frames
        % average frames for scan
        img_sum = 0;
        % THIS WILL FAIL IF YOUR FRAMES ARE NOT INDEXED FROM 0000
        for f = 0:n_frames-1
            % Create filename for sample
            filename = strcat(dir_2d, file_hdr, num2str(scan), sprintf('_%04d', f), ext);

            % open the image & append extra rows
            img = vertcat(double(imread(filename)), zeros(rows_needed, im_col));
            img = medfilt2(img, [3 3]);
            img_sum = img_sum + img;
        end
        img_list = {img_sum ./ n_frames};
    else
        % initialize a list to collect all the images for individual
        % analysis
        img_list = {};
        % loop through all frames for the current scan
        for f = 0:n_frames-1
            % Create filename for sample
            filename = strcat(dir_2d, file_hdr, num2str(scan), sprintf('_%04d', f), ext);

            % open the image & append extra rows
            img = vertcat(double(imread(filename)), zeros(rows_needed, im_col));
            img = medfilt2(img, [3 3]);
            % append image to list of images
            img_list = [img_list, img];
        end
    end
    
    % loop through processed images--this will just be the averaged image
    % if averaging is turned (i.e., avg_frames == true)
    for frame=0:length(img_list)-1
        img = img_list{frame+1};
        % subtract background
        if two_para_subt
            %_________________perform 2 para subtraction___________________________
            [alpha(k), beta(k)] = Two_Para_Subtraction(img, bkg_avg, q_low_min,...
                q_low_max, q_high_min, q_high_max, q_matrix);
            bkg_2param = bkg_avg;
            bkg_2param(bkg_2param>0) = alpha(k)*bkg_2param(bkg_2param>0) + beta(k);
            % subtract background
            img_subt = img - bkg_2param; 
        else
            % subtracted average of background scans of glass mold
            img_subt = img - bkg_avg;
        end

        % plot the  image
        if open_new_windows
            h = figure(k);  
        else
            clf
            h = figure(1);
        end
        % set range of colorbar
        imagesc (img_subt, [0, max_intensity])  
        axis image
        axis off
        colorbar
        colormap jet
        title_str = 'Bkgd-subt Intensity';
        if avg_frames
            title_str = sprintf('%s Scan %d (averaged)', title_str, scan);
        else
            title_str = sprintf('%s Scan %d Frame %d', title_str, scan, frame);
        end
        title(title_str)
        % add legend label for current sample
        if avg_frames
            labels = [labels, {sprintf('Scan %d', scan)}];
        else
            % list frame number with scan
            labels = [labels, {sprintf('Scan %d.%d', scan, frame)}];
        end
        % save the image
        if avg_frames
            frame_tag = '_avg';
        else
            frame_tag = sprintf('%04d', frame);
        end
        filename = strcat(dir_2d_bkgd_subt, 'S', num2str(scan),...
            '_', detector_type, frame_tag); 
        print (h, filename, '-dpng', '-r300'); 

        % 1D INTENSITIES
        %interpolate values for integrating in R or Phi for 1D intensities
        unwrapped_subt = interp2(x, y, img_subt, xi, yi);
        % clean up interpolation by removing nan and inf
        unwrapped_subt(isnan(unwrapped_subt)) = 0;
        unwrapped_subt(isinf(unwrapped_subt)) = 0;
        %calculate I(q) and I(phi)
        [q_1D, IvsQ_Sub(:,ct)] = IvsQ_Calculator(unwrapped_subt , row_center, col_center,...
            r_min, r_max, parameters);
        [phi_1D, IvsPhi_Sub(:,ct)] = IvsPhi_Calculator(unwrapped_subt, row_center, col_center,...
            r_min, r_max, parameters);

        % Repeat above process for raw images (background not subtracted)
        unwrapped = interp2(x, y, img, xi, yi);
        unwrapped(isnan(unwrapped)) = 0;
        unwrapped(isinf(unwrapped)) = 0;
        %calculate I(q) and I(phi)
        [~, IvsQ(:,ct+1)] = IvsQ_Calculator(unwrapped , row_center, col_center,...
            r_min, r_max, parameters);
        [~, IvsPhi(:,ct+1)] = IvsPhi_Calculator(unwrapped, row_center, col_center,...
            r_min, r_max, parameters);
        % update header, scan is ones place, frame is thousandths
        hdr(ct) = scan + 0.001*frame;
        % increment the counter
        ct = ct + 1;
    end
end

% add raw background to 1D intensities
unwrapped_bkg = interp2(x, y, bkg_avg, xi, yi);
unwrapped_bkg(isnan(unwrapped_bkg)) = 0;
unwrapped_bkg(isinf(unwrapped_bkg)) = 0;
% Save the background intensities
[~, IvsQ(:,1)] = IvsQ_Calculator(unwrapped_bkg , row_center, col_center,...
    r_min, r_max, parameters);
[~, IvsPhi(:,1)] = IvsPhi_Calculator(unwrapped_bkg, row_center, col_center,...
    r_min, r_max, parameters);

%% save data files of 1D intensities
% -1 signifies background
csvwrite (strcat(dir_data, 'q_1D.csv'), q_1D);
csvwrite (sprintf('%s%s_Sub_IvsQ.csv', dir_data, detector_type), vertcat(hdr, IvsQ_Sub));
csvwrite (sprintf('%s%s_IvsQ.csv', dir_data, detector_type), vertcat([-1, hdr], IvsQ));
csvwrite (strcat(dir_data, 'phi_1D.csv'), phi_1D);
csvwrite (sprintf('%s%s_Sub_IvsPhi.csv', dir_data, detector_type), vertcat(hdr, IvsPhi_Sub));
csvwrite (sprintf('%s%s_IvsPhi.csv', dir_data, detector_type), vertcat([-1, hdr], IvsPhi));

%% I vs. q (1D): plot averaged background next to corresponding raw data

% ignore nan values in background
not_nan = ~isnan(IvsQ(:,1));
% set axis limits
xmin = min(q_1D(not_nan));
xmax = max(q_1D(not_nan));
ymin = -40;
ymax = max(IvsQ(:));

if save_I_1D
    ct = 1;
    for k = 1:length(scan_arr)
        scan = scan_arr(k);
        scan_str = sprintf('Scan %d', scan);
        for frame=0:n_frames_2_save-1
            if open_new_windows
                figure(3500+ct)
            else
                figure(3500)
            end
            hold on
            plot (q_1D, IvsQ(:,ct+1), 'r-', 'linewidth', 2)
            if two_para_subt
                curr_bkg = alpha(ct)*IvsQ(:,1)+beta(ct);
                bkg_str = ' with 2-param-corrected bkg';
                bkg_label = 'corrected bkg';
                filename = sprintf('%sIvsQ_2param_bkg_scan_%d',dir_1d, scan);
            else
                curr_bkg = IvsQ(:,1);
                bkg_str = ' with averaged bkg';
                bkg_label = 'avg bkg';
                filename = sprintf('%sIvsQ_avg_bkg_scan_%d',dir_1d, scan);
            end
            if avg_frames
                frame_str = ' (avg)';
            else
                frame_str = sprintf(' Frame %d', frame);
            end
            plot(q_1D, curr_bkg, 'b--', 'linewidth', 2)
            axis ([xmin xmax ymin ymax])
            xlabel ('q(1/Å)')
            ylabel ('Intensity[a.u.]')
            title_str = strcat(scan_str, frame_str, bkg_str);
            title(title_str)
            legend('raw data', bkg_label)
            set(legend,'location','best')
            set(gcf,'PaperPositionMode','auto');
            box on
            h=findall(gcf,'tag','legend'); % find the handle to the legend
            print(gcf, filename, '-dpng', '-r300'); 
            % increment counter
            ct = ct + 1;
        end
    end
end

%% I vs. phi (1D): plot averaged background next to corresponding raw data

% set axis limits
xmin = phi_min;
xmax = phi_max;
ymin = -40;
ymax = max(IvsPhi(:));

if save_I_1D
    ct = 1;
    for k = 1:length(scan_arr)
        for frame=0:n_frames_2_save-1
            scan = scan_arr(k);
            scan_str = sprintf('Scan %d', scan);
            if open_new_windows
                figure(4500+ct)
            else
                figure(4500)
            end
            hold on
            plot (phi_1D, IvsPhi(:,ct+1), 'r-', 'linewidth', 2)
            if two_para_subt
                curr_bkg = alpha(ct)*IvsPhi(:,1)+beta(ct);
                bkg_str = ' with 2-param-corrected bkg';
                bkg_label = 'corrected bkg';
                filename = sprintf('%sIvsPhi_2param_bkg_scan_%d',dir_1d, scan);
            else
                curr_bkg = IvsPhi(:,1);
                bkg_str = ' with averaged bkg';
                bkg_label = 'avg bkg';
                filename = sprintf('%sIvsPhi_avg_bkg_scan_%d',dir_1d, scan);
            end
            if avg_frames
                frame_str = ' (avg)';
            else
                frame_str = sprintf(' Frame %d', frame);
            end
            title_str = strcat(scan_str, frame_str, bkg_str);
            title(title_str)
            plot(phi_1D, curr_bkg, 'b--', 'linewidth', 2)
            axis ([xmin xmax ymin ymax])
            xlabel ('phi [deg]')
            ylabel ('Intensity[a.u.]')
            title(title_str)
            legend('raw data', bkg_label)
            set(legend,'location','best')
            set(gcf,'PaperPositionMode','auto');
            box on
            h=findall(gcf,'tag','legend'); % find the handle to the legend
            print(gcf, filename, '-dpng', '-r300');
            % increment counter
            ct = ct + 1;
        end
    end
end

%% plot the change in alpha & beta 
if two_para_subt
    figure(4003)
    plot (alpha, '--ok')
    xlabel ('Frame'); 
    ylabel ('Alpha value')
    title('Alpha Parameter vs. Frame')

    figure(4004)
    plot (beta, '--ok')
    xlabel ('Frame'); 
    ylabel ('Beta value')
    title('Beta Parameter vs. Frame')
end
