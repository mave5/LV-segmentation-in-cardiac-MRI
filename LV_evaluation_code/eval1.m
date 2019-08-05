%eval.m
%   Demonstration of how to use compare_contours.m and write_results.m and
%   evaluate a batch of contours.
%
%   Copyright: Imaging Research, Sunnybrook Health Sciences Centre, Toronto, ON, Canada
%   Author: Perry Radau and Yingli Lu
%   Email: perry.radau@gmail.com
%   Version: 1.0
%   Date: 2009/06/29
clc
clear all
warning off all
addpath 'functions'
%all studies with SC* prefix under the root_folder directory will be processed.
folder_prefix = 'SC*'; %This prefix is used for all studies.
process_all = true; %Set this to false if you want to only calc results for first two studies.
%DICOM folder name
para.dicom_foldername = 'DICOM';
%manual contour folder name
para.manual_contour_foldername = 'contours-manual';
para.manual_contour_subfoldername = 'IRCCI-expert';
%auto contour folder name
para.auto_contour_foldername = 'contours-auto';
para.auto_contour_subfoldername = 'auto2';
%para.auto_contour_subfoldername = 'BIL_HS';

%If set to 1,normal is calculated based on auto contours, manual constour is the reference
%contour. If set to 0, normal is calculated based on manual contours, auto contour is the reference contour
para.auto_based_noraml = 1; 

%distance figure folder name
para.distance_figure_foldername = 'distance_figure';
%DIDOM file name prefix 
para.name_prefix = 'IM-0001-'; %according to the IM-0001-XXXX naming standard
para.digit_length=4; %IM-0001-XXXX, XXXX is 4 digits
%contour mode
para.inside_contour_mode = 'icontour';
para.outside_contour_mode = 'ocontour';
%save distance figure if set to 1
para.save_figure = 1;
%segmentation mode
para.auto_seg_mode='auto';
para.manual_seg_mode ='manual';
%for average distance
para.init_value = -999;
para.dist_limit = 5; %mm. If a contour's average distance larger than this number, it will be considered as bad contour

%Keep checking directory until valid
count = 0;
while (true)
    %root_folder = input('Enter root folder: ', 's');
    root_folder=uigetdir;
    results_folder = fullfile(root_folder,'results');
    dirs = dir([root_folder filesep folder_prefix]);
    if ~isempty(dirs)
        break;
    end
    %check path for input files
    disp('Invalid root folder, retry.');
    count = count + 1;
    if (count >= 3)
        return;
    end
end

%initialize
all_compare_result = cell(length(dirs),1);

if (process_all)
    num_studies = length(dirs);
else 
    num_studies =3;%use this for testing.
end
    
for i = 1:num_studies
    display(['Processing: ' dirs(i).name]);
    %full dicom folder. The DICOM images must be in this folder.
    dicom_folder = [root_folder filesep dirs(i).name filesep para.dicom_foldername];
    % 
    %full manual contour folder.The last subfolder indicates which expert
    %drew the manual contours.
    manual_contour_folder = [root_folder filesep dirs(i).name filesep para.manual_contour_foldername filesep para.manual_contour_subfoldername];
    
    if ~exist(manual_contour_folder,'dir')
           display(['Missing manual contour folder:', manual_contour_folder]);
           return;
    end
    
    %full auto contour folder.The last subfolder indicates which algorithm
    %was used.
    auto_contour_folder = [root_folder filesep dirs(i).name filesep para.auto_contour_foldername filesep para.auto_contour_subfoldername];
    if ~exist(auto_contour_folder,'dir')
           display(['Missing auto contour folder:', auto_contour_folder]);
           continue;
    end
    %compare
    compare_result = compare_contours(dicom_folder, manual_contour_folder, auto_contour_folder,para);
    %record results. 
    %The '_i' indicates inner contour, 
    %The '_o' indicates outer contour. 
    %The '_pic' means Papillary Included in the LV Cavity,
    %The '_pim' means Papillary Included in the Myocardium

    if ~isempty(compare_result)
        all_compare_result{i}.patient_id = dirs(i).name;
        all_compare_result{i}.auto_number_i = compare_result.auto_number_i;
        all_compare_result{i}.auto_number_o=  compare_result.auto_number_o;
        all_compare_result{i}.manual_number_i = compare_result.manual_number_i;
        all_compare_result{i}.manual_number_o = compare_result.manual_number_o;
        all_compare_result{i}.detect_percent_i = compare_result.detect_percent_i;
        all_compare_result{i}.detect_percent_o = compare_result.detect_percent_o;
        all_compare_result{i}.good_percent_i = compare_result.good_percent_i;
        all_compare_result{i}.good_percent_o = compare_result.good_percent_o;
        
        all_compare_result{i}.auto_ef_pic = compare_result.auto_ef_pic;
        all_compare_result{i}.auto_ef_pim = compare_result.auto_ef_pim;
        all_compare_result{i}.auto_lvm_pic= compare_result.auto_lvm_pic;
        all_compare_result{i}.auto_lvm_pim= compare_result.auto_lvm_pim;
        
        all_compare_result{i}.manual_lvm_pic = compare_result.manual_lvm_pic;
        all_compare_result{i}.manual_lvm_pim = compare_result.manual_lvm_pim;
        all_compare_result{i}.manual_ef_pic = compare_result.manual_ef_pic;
        all_compare_result{i}.manual_ef_pim = compare_result.manual_ef_pim;
        
        all_compare_result{i}.avg_dist_i = compare_result.avg_dist_i;
        all_compare_result{i}.avg_dist_o = compare_result.avg_dist_o;
        all_compare_result{i}.avg_dm_i = compare_result.avg_dm_i;
        all_compare_result{i}.avg_dm_o = compare_result.avg_dm_o;

        % manual and auto volumes
        all_compare_result{i}.auto_edv_pic= compare_result.auto_edv_pic;
        all_compare_result{i}.auto_edv_pim= compare_result.auto_edv_pim;       
        all_compare_result{i}.manual_edv_pic = compare_result.manual_edv_pic;
        all_compare_result{i}.manual_edv_pim = compare_result.manual_edv_pim;

        all_compare_result{i}.auto_esv_pic= compare_result.auto_esv_pic;
        all_compare_result{i}.auto_esv_pim= compare_result.auto_esv_pim;       
        all_compare_result{i}.manual_esv_pic = compare_result.manual_esv_pic;
        all_compare_result{i}.manual_esv_pim = compare_result.manual_esv_pim;
       
   end
end

%write result to xls(excel) file
write_results(results_folder, all_compare_result);

disp('Done!');
