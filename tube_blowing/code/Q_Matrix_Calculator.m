function [Q_2D] = Q_Matrix_Calculator()

% the "truncated" image size is 981x981
im_row = 978; 
im_col = 1920; 

row_center = 978;
col_center = 958;

%beamtime parameters
lamda = 0.7293*10^-10; % in m
pixsize = 88.6 *10^-6; % in m  
SaDet = 199.87 * 10^-3 ; % in m

%create the radius matrix
for i = 1:im_row
    for j = 1:im_col
        row_distance = abs(i - row_center);
        col_distance = abs(j - col_center);
        radius(i,j) = sqrt((row_distance^2) + (col_distance^2));
    end
end

%create a q matriz
Q_2D = 4 * 10^-10*pi*sin(0.5*atan(radius*pixsize/SaDet))/lamda; 

% plot
% figure(5000)
% imagesc (Q_2D)
% axis off
% axis image
% colormap jet
% colorbar

end

