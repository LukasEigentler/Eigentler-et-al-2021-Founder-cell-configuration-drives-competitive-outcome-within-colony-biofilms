function [blue_length,accept] = AFS(filepath, identifier,size_threshold)
% AFS  Computes access to free space score based on images
%   [blue_length,accept] = AFS(filepath, identifier,size_threshold)
%   computes the access to free space score blue_length based on the images
%   containing identifier in their filenames. Coloured regiions of size
%   smaller than size_threshold are rejected as noise. The function also
%   returns accept, a user-defined string that enables the parent code to
%   either accept or reject the outcome.
%
% Author: Lukas Eigentler
% Last updated: 08/06/2021

close all;
Files = dir(filepath); % read all files
files_ind = find(contains({Files.name},identifier+"]")==1); % find files containing identifier
blue_image_ind = find(contains({Files(files_ind).name},'blue')==1); % find which of remaining files contain blue
green_image_ind = find(contains({Files(files_ind).name},'green')==1); % find which of remaining files contain green
merged_image_ind = find(contains({Files(files_ind).name},'merged')==1); % find which of remaining files contain merged
R = 1000; % radius of reference circle (in pixels)
%% determine coloured patches
[~, binaryImage_blue] = DeltaE_mod([filepath,Files(files_ind(blue_image_ind)).name],10); % determine matching colours blue
 sz = size(binaryImage_blue);
CC_blue = bwconncomp(binaryImage_blue); % find connected blue regions
region_size_blue = cellfun(@numel,CC_blue.PixelIdxList); % find sizes of regions
reject_ind_blue = find(region_size_blue < size_threshold); % reject regions that are too small as noise
accept_ind_blue = find(region_size_blue >= size_threshold); % accept regions that are big enough
 for rr = 1:length(reject_ind_blue)
    binaryImage_blue(CC_blue.PixelIdxList{reject_ind_blue(rr)}) = 0; % remove rejected regions from binary image
 end
 blue_ind = [0,0];
 for aa = 1:length(accept_ind_blue) % get indices of all pixels in accepted regions
     [blue_ind_y,blue_ind_x] = ind2sub(sz,CC_blue.PixelIdxList{accept_ind_blue(aa)});
     blue_ind_new = [blue_ind_x, blue_ind_y];
     blue_ind = [blue_ind; blue_ind_new];
 end
 blue_ind = blue_ind(2:end,:);
 
[~, binaryImage_green] = DeltaE_mod([filepath,Files(files_ind(green_image_ind)).name],40); % determine matching colours green
CC_green = bwconncomp(binaryImage_green); % find connected green regions
region_size_green = cellfun(@numel,CC_green.PixelIdxList); % find sizes of regions
reject_ind_green = find(region_size_green < size_threshold); % reject regions that are too small as noise
accept_ind_green = find(region_size_green >= size_threshold); % accept regions that are big enough
 for rr = 1:length(reject_ind_green)
    binaryImage_green(CC_green.PixelIdxList{reject_ind_green(rr)}) = 0; % remove rejected regions from binary image
 end
 
green_ind = [0,0];
 for aa = 1:length(accept_ind_green) % get indices of all pixels in accepted regions
     [green_ind_y,green_ind_x] = ind2sub(sz,CC_green.PixelIdxList{accept_ind_green(aa)});
     green_ind_new = [green_ind_x, green_ind_y];
     green_ind = [green_ind; green_ind_new];
 end
 green_ind = green_ind(2:end,:);
 
 binaryImage_combined = or(binaryImage_blue,binaryImage_green); % create combined binary image
 
 %% find centre

 x = 1:sz(2); y = 1:sz(1);
 total_mass = sum(binaryImage_combined(:)); % calc total mass
 [row,col] = find(binaryImage_combined == 1); %find loc of mass
 centre_of_mass = [0,0];
 for rr = 1:length(row) % find centre of mass
     centre_of_mass = centre_of_mass + [row(rr),col(rr)];
 end
 centre_of_mass = centre_of_mass/total_mass;
blue_ind = blue_ind(vecnorm(blue_ind-fliplr(centre_of_mass),2,2)<R/2,:);
green_ind = green_ind(vecnorm(green_ind-fliplr(centre_of_mass),2,2)<R/2,:);

% recalc centre after excluding outside noise
col_new = [blue_ind(:,1);green_ind(:,1)];
row_new = [blue_ind(:,2);green_ind(:,2)];
total_mass_new = length(row_new); % calc total mass
centre_of_mass_new = [0,0];
for rr = 1:length(row_new) % find centre of mass
    centre_of_mass_new = centre_of_mass_new + [row_new(rr),col_new(rr)];
end
centre_of_mass = centre_of_mass_new/total_mass_new;

 %% plotting
 figure

col = lines;
green = [0,1,0];
mag = [0, 1, 1];
theta_col = linspace(0,2*pi,1000);

subplot(1,2,1)
plot(blue_ind(:,1), blue_ind(:,2),'.', 'color', mag)
axis equal
axis on
hold on
plot(green_ind(:,1), green_ind(:,2),'.', 'color', green)
xlim([min(centre_of_mass(2) + R*cos(theta_col)),max(centre_of_mass(2) + R*cos(theta_col))])
ylim([min(centre_of_mass(1) + R*sin(theta_col)),max(centre_of_mass(1) + R*sin(theta_col))])


subplot(1,2,2)
imshow(imread([filepath,Files(files_ind(merged_image_ind)).name]))
axis on
hold on
xlim([min(centre_of_mass(2) + R*cos(theta_col)),max(centre_of_mass(2) + R*cos(theta_col))])
ylim([min(centre_of_mass(1) + R*sin(theta_col)),max(centre_of_mass(1) + R*sin(theta_col))])

dtheta = 2*pi/1000;
blue_length = 0;
for theta = 0:dtheta:2*pi-dtheta
    [~,dist_blue] = dsearchn(blue_ind,[centre_of_mass(2) + R*cos(theta), centre_of_mass(1) + R*sin(theta)]); %calc distance to clostest blue patch
    [~,dist_green] = dsearchn(green_ind,[centre_of_mass(2) + R*cos(theta), centre_of_mass(1) + R*sin(theta)]); %calc distance to clostest green patch
    [~,min_ind] = min([double(dist_blue),double(dist_green)]); % calc which colour patch is closer
    if min_ind == 1
        subplot(1,2,1)
        plot(centre_of_mass(2) + R*cos(theta), centre_of_mass(1) + R*sin(theta), '.','MarkerSize',10, 'color', mag)
        subplot(1,2,2)
        plot(centre_of_mass(2) + R*cos(theta), centre_of_mass(1) + R*sin(theta), '.','MarkerSize',10, 'color', mag)
        blue_length = blue_length + dtheta;
    else
        subplot(1,2,1)
        plot(centre_of_mass(2) + R*cos(theta), centre_of_mass(1) + R*sin(theta), '.','MarkerSize',10, 'color', green)
        subplot(1,2,2)
        plot(centre_of_mass(2) + R*cos(theta), centre_of_mass(1) + R*sin(theta), '.','MarkerSize',10, 'color', green)
    end
end
subplot(1,2,1)
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
accept = input('Accept or reject and continue or reject outcome? (a/c/r)', 's');
if accept == 'a'
    saveas(gcf, ['figures/',strrep(strrep(strrep(Files(files_ind(merged_image_ind)).name, '.lif','_'), '.tiff', '_matlab_figure'), '.jpg', '_matlab_figure')])
end
%%
