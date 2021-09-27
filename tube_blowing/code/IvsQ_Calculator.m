function [q_1D, IvsQ] = IvsQ_Calculator (image, bc_row, bc_col, r_min,...
    r_max, parameters)

    % This script generates IvsQ and IvsPhi (phi = azimuthal angle on image) by
    % "unwrapping" an image, that is transforming the pattern from polar
    % coordinates to cartesian coordinates by interpolation on a grid of chosen
    % dimensions. The summation along the x axis
    % gives the IvsQ and the summation along the theta axis gives IvsPhi

    set(0,'defaultAxesFontSize',18);
    set(0,'defaultTextFontSize',18);
    set(0,'defaultTextFontName','Times');
    set(0,'defaultTextFontWeight','bold');
    
    % the beamtime parameters
    lamda = parameters(1,1);
    pixsize = parameters(2,1);
    SaDet = parameters(3,1); 
    Np = size(image,1);
    xc = bc_col; % corrected values from BeamCenter code
    yc = bc_row;  % corrected values from BeamCenter code
    w = Np * pixsize;

    % the parameters that set the number of points and the range of the summation
    nt = round(2*pi*r_max); % number of theta points
    nr = r_max-r_min;
    th1 = -pi/2;
    th2 = 3*pi/2; 
    
    %read in the image file
    im = image;
    
    %create the mesh in x, y 
    [x,y] = meshgrid(1:size(im,2),1:size(im,1)); 

    % create the mesh in phi and r
    [phi, radius] = meshgrid(linspace(th1,th2,nt), linspace(r_min,r_max,nr));
%     xi = xc+radius.*cos(phi); % these are the cartesian coordinates of each pixel
%     yi = yc+radius.*sin(phi);
% 
%     % 2D interpolation: interpolate the values at (xi, yi) of the original
%     % the image im so that each pixel corresponds to an (xi, yi) of the unwrapped image
%     mapping = interp2(x, y, double(im), xi, yi); 
%     disp(max(max(mapping)))
%     mapping(isnan(mapping)) = 0; % fill out the NaN values with zero
% 
%     disp(max(max(mapping)))
    % The averages give the IvsQ and IvsPhi
    IvsQ_tot = sum(im, 2);
    num_live_pixels = sum(im > 0, 2);
    IvsQ = IvsQ_tot./num_live_pixels;
    
%     % plot the unwrapped image
%     figure(1000)
%     imagesc (mapping, [0 100])
%     colormap jet
%     colorbar

    % obtaining the q_matrix
    qq = 4 * 10^-10*pi*sin(0.5*atan(radius*pixsize/SaDet))/lamda; % the Q grid
    q_1D = qq(:,1); %the q range for this IvsQ, confirm that it's the same length for all analysis
    phi = phi(1,:)*180/pi; % azimuthal angle on detector

end