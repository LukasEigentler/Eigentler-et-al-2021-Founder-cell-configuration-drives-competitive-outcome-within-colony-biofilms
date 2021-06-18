%% Access to free space score image analysis - non-isogenic pair

% This code calculates the access to free space score from an image for the
% non-isogenic strain pair. For the corresponding code for the isogenic
% strain pair, see "image_analysis_theta_score.m".
% For each case to be analysed, the script requires 3 image files (.jpg,
% .tiff, etc) in the folder "unprocessed": two single-channel images and 
% one image showing both channels. All three files require to have a unique
% identifier within a pair of square brackets in their filename. Moreover,
% the filenames of the single-channel images require the words "blue" and 
% "green", respectively, and the merged image requires the word "Merged" to
% appear in its filename.
% The script is run interactively. The user is presented with both
% single-strain images and is required to select a coloured area in each
% image. The algorithm then selects all similarly coloured regions using
% Delta E colour difference, calculates the access to free space score,
% and displays a plot of identified colour regions next to the merged
% image. The user is required to check the result for accuarcy. If
% approved, the file writes data to a file and moves the processed image
% files into the folder "Processed". The constructed Matlab plot alongside
% the "Merged" image is saved in the folder "Figures".
% Please note the license file for the function "DeltaE_mod.m".

% Author: Lukas Eigentler
% Last updated: 08/06/2021

clear; close all;
% cas = '3610_vs_6153'; % select strain pair
% cas = '6153';
cas = '3610';
size_threshold = 500; % minimum size of element to be detected (smaller elements will be rejected as noise)
unprocessed = dir('unprocessed'); % read image files to be processed
identifier = strings(length(unprocessed)-2,1); % find unique identifiers of images
for file_ind = 3:length(unprocessed)
    identifier_start = strfind(unprocessed(file_ind).name, '[');
    identifier_end = strfind(unprocessed(file_ind).name, ']');
    identifier(file_ind-2) = unprocessed(file_ind).name(identifier_start+1:identifier_end-1);
end
identifier = unique(identifier);

%% load existing data/check if data exists
if isfile(['initial_theta_score_image_analysis_results_', cas,'.csv'])
    results = readtable(['initial_theta_score_image_analysis_results_', cas,'.csv'], 'Delimiter',',');
else
    Filename = "start"; blue_theta_length = NaN;
    results = table(Filename, blue_theta_length);
end


%% process all images
for file_ind = 1:min([100,length(identifier)])
    
    accept = 'c';
    while accept == 'c'
        fprintf("Analysing "+identifier(file_ind))
        [blue_length, accept] = AFS('unprocessed/',identifier(file_ind),size_threshold); % image processing
        if accept == 'a'
            move_ind = find(contains({unprocessed.name},identifier(file_ind)+"]")==1); % find files containing identifier
            end_ind = strfind(unprocessed(move_ind(1)).name,']'); % determine filename to be written to results file 
            fname = unprocessed(move_ind(1)).name;
            fname = fname(1:end_ind);
            for ff = 1:length(move_ind)
                movefile(['unprocessed/',unprocessed(move_ind(ff)).name],['processed/',unprocessed(move_ind(ff)).name]) % move files to processed folder
            end  
         
            if isnan(results.blue_theta_length(1))
                l = 0;
                new_row = table(string(fname), blue_length);
                results(l+1,:) = new_row; % add new results to table and delete placeholder
            else
                
                new_row = table(cellstr(string(fname)), blue_length);
                if isempty(find(strcmp(string(fname), results.Filename) == 1, 1)) % check if image is already in results
                    l = length(results.Filename);
                    results(l+1,:) = new_row; % add new results to existing table
                else
                    l = find(strcmp(string(fname), results.Filename) == 1);
                    results(l,:) = new_row; % add new results to existing table
                end
            end

        elseif accept == 'c'
            size_threshold = input(['Enter new size threshold. Current value is', num2str(size_threshold)]);
        end
    end
    
end
writetable(results, ['initial_theta_score_image_analysis_results_', cas,'.csv']) % save results

