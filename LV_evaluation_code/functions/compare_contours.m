function compare_result = compare_contours(dicom_path,manual_contour_path,auto_contour_path,para)
%COMPARE_CONTOURS Compare manual rawn contours with auto contours
%   COMPARE_CONTOURS(DICOM_PATH,MANUAL_CONTOUR_PATH,AUTO_CONTOUR_PATH)
%
%   Copyright: Imaging Research, Sunnybrook Health Sciences Centre, Toronto, ON, Canada
%   Author: Perry Radau and Yingli Lu
%   Email: perry.radau@gmail.com
%   Version: 1.0
%   Date: 2009/06/29

%   OUTPUT: compare_result is a struct with fields of
%-----------------------------------------------------
%   auto_number_i:  total number of auto inside contours
%   auto_number_o:  total number auto outside contours;
%   manual_number_i: total number manual inside contours
%   manual_number_o: total number manual outside contours
%   detect_percent_i: percent of detected number of auto inside contours compared to total number of  manual inside contours
%   detect_percent_o: percent of detected number of auto outside contours compared to total number of  manual outside contours
%   good_percent_i: percent of good auto inside contours (a good contour has an average distance smaller than para.dist_limit)
%   good_percent_o: percent of good auto outside contours
%   auto_ef_pic: auto contour's ejection fraction, '_pic' means Papillary Included in the LV Cavity
%   auto_ef_pim: auto contour's ejection fraction, '_pim' means Papillary Included in the Myocardium
%   auto_lvm_pic: auto contour's lv mass
%   auto_lvm_pim: 
%   manual_ef_pic: auto contour's ejection fraction
%   manual_ef_pim: 
%   manual_lvm_pic: manual contour's lv mass
%   manual_lvm_pim: 
%   avg_dist_i: average distance of inside contours
%   avg_dist_o: average distance of outside contours
%   avg_dm_i: average dice metric of inside contours
%   avg_dm_o: average dice metric of outside contours

%   auto_missing_index_i: missing auto inside contours' index
%   auto_missing_index_o: missing auto outside contours' index
%   auto_bad_index_i: bad auto inside contours' index (a bad contour has an  average distance larger than para.dist_limit)
%   auto_bad_index_o: bad auto outside contours' index 
%-----------------------------------------------------

%-initialize
compare_result = [];

%-check path
if ~exist(dicom_path,'dir')
    disp([dicom_path ' :NOT exist!'])
    return;
end
if ~exist(manual_contour_path,'dir')
    disp([manual_contour_path ' :NOT exist!'])
    return;
end
if ~exist(auto_contour_path,'dir')
    disp([auto_contour_path ' :NOT exist!'])
    return;
end

%-check files
dicom_files = dir([dicom_path filesep '*.dcm']);
if isempty(dicom_files)
    disp([dicom_path ' :NO dicom files!'])
    return;
end
manual_contour_files = dir([manual_contour_path filesep '*.txt']);
if isempty(manual_contour_files)
    disp([manual_contour_path ' :NO manual contours!'])
    return;
end
auto_contour_files = dir([auto_contour_path filesep '*.txt']);
if isempty(auto_contour_files)
    disp([auto_contour_path ' :NO auto contours!'])
    return;
end

%-dicominfo
try
  dicom_filename = dicom_files(1).name; %use the first dicom file.
  full_dicom_filename = [dicom_path filesep dicom_filename];
  dicom_meta= dicominfo(full_dicom_filename);
  
  para.width = dicom_meta.Width;%image width
  para.height = dicom_meta.Height;%image height
  para.pixel_spacing = dicom_meta.PixelSpacing; % mm
  para.thickness = dicom_meta.SliceThickness;% mm
  para.gap = dicom_meta.SpacingBetweenSlices - para.thickness; % mm
  para.phase_number = dicom_meta.CardiacNumberOfImages; %numer of phases
  para.image_number = length(dicom_files);%number of images
catch
  s = lasterror;
  disp(s.message);
  return
end

%-calc percent of auto contours compared to manual contours
manual_number_i = 0; %number of manual inside contours
manual_number_o = 0; %number of manual outside contours
manual_index_i = []; %index of manual inside contours
manual_index_o = []; %index of manual outside contours

%count manual contours
for ix = 1:para.image_number 
    sindex = add_zero_index(ix,para.digit_length);
    %inside contour
    full_icontour_filename = get_contour_filename(manual_contour_path, para.name_prefix, sindex,para.inside_contour_mode,para.manual_seg_mode);
    if exist(full_icontour_filename,'file') 
       fileinfo_i = dir(full_icontour_filename);
       if fileinfo_i.bytes ~= 0 % not null
          manual_number_i = manual_number_i + 1;
          manual_index_i = [manual_index_i; ix]; %#ok<AGROW>
       end
    end
    %outside contour
    full_ocontour_filename = get_contour_filename(manual_contour_path,para.name_prefix, sindex, para.outside_contour_mode,para.manual_seg_mode);
    if exist(full_ocontour_filename,'file') 
       fileinfo_o = dir(full_ocontour_filename);
       if fileinfo_o.bytes ~= 0 % not null
          manual_number_o = manual_number_o + 1;
          manual_index_o = [manual_index_o; ix]; %#ok<AGROW>
       end
    end
end

if (manual_number_i == 0)
   disp([manual_contour_path ' :NO inside contours!'])
   return; 
end
if (manual_number_o == 0)
   disp([manual_contour_path ' :NO outside contours!'])
   return; 
end

para.manual_index_i = manual_index_i;
para.manual_index_o = manual_index_o;

%count auto inside contours according to manual_index_i 
auto_number_i = 0;
auto_missing_index_i = [];
for ix = 1:length(manual_index_i)
    sindex = add_zero_index(manual_index_i(ix),para.digit_length);
    full_icontour_filename = get_contour_filename(auto_contour_path,para.name_prefix, sindex, para.inside_contour_mode, para.auto_seg_mode);
    if exist(full_icontour_filename,'file') 
       fileinfo_i = dir(full_icontour_filename);
       if fileinfo_i.bytes ~= 0 % not null
          auto_number_i = auto_number_i + 1;  
       else
          auto_missing_index_i = [auto_missing_index_i  manual_index_i(ix)]; %#ok<AGROW> 
       end  
    else
       auto_missing_index_i = [auto_missing_index_i  manual_index_i(ix)]; %#ok<AGROW>
    end
end

if (auto_number_i == 0)
   disp([auto_contour_path  ' :NO inside contours!'])
   return; 
end

%count auto outside contours according to manual_index_o 
auto_number_o = 0;
auto_missing_index_o = [];
for ix = 1:length(manual_index_o)
    sindex = add_zero_index(manual_index_o(ix),para.digit_length);
    full_ocontour_filename = get_contour_filename(auto_contour_path,para.name_prefix, sindex, para.outside_contour_mode,para.auto_seg_mode);
    if exist(full_ocontour_filename,'file') 
       fileinfo_o = dir(full_ocontour_filename);
       if fileinfo_o.bytes ~= 0 % not null
          auto_number_o = auto_number_o + 1;
       else
           auto_missing_index_o = [auto_missing_index_o  manual_index_o(ix)]; %#ok<AGROW>
       end
    else
       auto_missing_index_o = [auto_missing_index_o  manual_index_o(ix)]; %#ok<AGROW>
    end
end

if (auto_number_o == 0)
   disp([auto_contour_path ' :NO outside contours!'])
   %return; 
end

%-calc perpendicular distance & Dice Metric
%inside contours
avg_dist_i=ones(length(manual_index_i),1)*para.init_value;
dm_i = ones(length(manual_index_i),1)*para.init_value;
auto_bad_index_i = [];

for ix = 1:length(manual_index_i)
    sindex = add_zero_index(manual_index_i(ix),para.digit_length);
    % manual contour
    manual_icontour_filename = get_contour_filename(manual_contour_path,para.name_prefix,sindex, para.inside_contour_mode,para.manual_seg_mode);
    % auto contour
    auto_icontour_filename = get_contour_filename(auto_contour_path,para.name_prefix,sindex, para.inside_contour_mode,para.auto_seg_mode);
    % if both exist, compare
    if (exist(manual_icontour_filename,'file') &&  exist(auto_icontour_filename,'file'))
         manual_xy =  load(manual_icontour_filename);
         auto_xy =  load(auto_icontour_filename);
         
         max_auto_x = max(auto_xy(:,1));
         max_auto_y = max(auto_xy(:,2));
         min_auto_xy = min(auto_xy(:));
         
         if (size(auto_xy,1) > 15  && min_auto_xy > 1 && max_auto_x < para.width && max_auto_y < para.height) %% if "point number of auto contour" < 1, may lead to difficulty to find a subset of reference that correspond to the current target
            try
                %perpendicular distance
                if (para.auto_based_noraml == 1) %calc normal based on auto contours, manual constour is the reference contour
                   avg_dist_i(ix) =calc_dist(auto_xy,manual_xy,auto_icontour_filename,para);
                else
                   avg_dist_i(ix) =calc_dist(manual_xy,auto_xy,auto_icontour_filename,para); 
                end
                
                %dice metric
                dm_i(ix) = calc_dm(auto_xy,manual_xy,para);
            catch
                s = lasterror;
                disp(s.message);
                continue;
            end
         end
     end
end

if max(avg_dist_i) == para.init_value
%    return;
end
auto_bad_ix_i =find(avg_dist_i >= para.dist_limit);
if ~isempty(auto_bad_ix_i)
    auto_bad_index_i = manual_index_i(auto_bad_ix_i);
end
auto_good_ix_i =find(avg_dist_i < para.dist_limit & avg_dist_i ~= para.init_value);
auto_good_index_i = manual_index_i(auto_good_ix_i);
para.auto_good_index_i = auto_good_index_i ;

%outside contours
auto_bad_index_o = [];
avg_dist_o=ones(length(manual_index_o),1)*para.init_value;
dm_o = ones(length(manual_index_o),1)*para.init_value;

for ix = 1:length(manual_index_o)
    sindex = add_zero_index(manual_index_o(ix),para.digit_length);
    % manual contour
    manual_ocontour_filename = get_contour_filename(manual_contour_path,para.name_prefix,sindex, para.outside_contour_mode,para.manual_seg_mode);
    % auto contour
    auto_ocontour_filename = get_contour_filename(auto_contour_path,para.name_prefix,sindex, para.outside_contour_mode,para.auto_seg_mode);
    % if both exist, compare
    if (exist(manual_ocontour_filename,'file') &&  exist(auto_ocontour_filename,'file'))
         manual_xy =  load(manual_ocontour_filename);
         auto_xy =  load(auto_ocontour_filename);
         
         max_auto_x = max(auto_xy(:,1));
         max_auto_y = max(auto_xy(:,2));
         min_auto_xy = min(auto_xy(:));
         
         if (size(auto_xy,1) > 15  && min_auto_xy > 1 && max_auto_x < para.width && max_auto_y < para.height)  % if "point number of auto contour" < 1, may lead to difficulty to find a subset of reference that correspond to the current target
            try
                %perpendicular distance
                if (para.auto_based_noraml == 1) %calc normal based on auto contours, manual constour is the reference contour
                   avg_dist_o(ix) =calc_dist(auto_xy,manual_xy,auto_ocontour_filename,para);
                else
                   avg_dist_o(ix) =calc_dist(manual_xy,auto_xy,auto_ocontour_filename,para); 
                end
                
                %dice metric
                dm_o(ix) = calc_dm(auto_xy,manual_xy,para);
            catch
               s = lasterror;
               disp(s.message);
               continue;
            end
         end
     end
end

if max(avg_dist_o) == para.init_value
    %return;
end
auto_bad_ix_o =find(avg_dist_o >= para.dist_limit);
if ~isempty(auto_bad_ix_o)
    auto_bad_index_o = manual_index_o(auto_bad_ix_o);
end
% if auto_good_ix_o is empty, mean(dm_o(auto_good_ix_o)) can be a NAN (undefined)
auto_good_ix_o = find(avg_dist_o < para.dist_limit & avg_dist_o ~= para.init_value);
auto_good_index_o = manual_index_o(auto_good_ix_o);
para.auto_good_index_o = auto_good_index_o;

%-calc ejection fraction and lv mass
para.auto = false;
lv_manual = calc_clinical_para(dicom_path, manual_contour_path, para);

para.auto = true;
lv_auto = calc_clinical_para(dicom_path, auto_contour_path, para,1,para.image_number/para.phase_number,lv_manual.es,lv_manual.ed);

%record result
compare_result.auto_number_i = auto_number_i; %total auto inside contours
compare_result.auto_number_o = auto_number_o;
compare_result.manual_number_i = manual_number_i;
compare_result.manual_number_o = manual_number_o;
compare_result.detect_percent_i = auto_number_i / manual_number_i * 100; %percent of detected auto inside contours
compare_result.detect_percent_o = auto_number_o / manual_number_o * 100;
compare_result.auto_missing_index_i = auto_missing_index_i; % missing auto inside contours
compare_result.auto_missing_index_o = auto_missing_index_o;

compare_result.auto_bad_index_i = auto_bad_index_i; %bad auto inside contours' index (average distance larger than para.dist_limit)
compare_result.auto_bad_index_o = auto_bad_index_o;
compare_result.good_percent_i = (auto_number_i - length(auto_bad_index_i)) /manual_number_i * 100; %percent of good auto inside contours
compare_result.good_percent_o = (auto_number_o - length(auto_bad_index_o)) /manual_number_o * 100;

compare_result.auto_ef_pic = lv_auto.ef_pic;%auto's ef, '_pic' means Papillary Included in the LV Cavity 
compare_result.auto_ef_pim = lv_auto.ef_pim; %'_pim' means Papillary Included in the Myocardium
compare_result.auto_lvm_pic = lv_auto.lvm_pic; %auto's lv mass
compare_result.auto_lvm_pim = lv_auto.lvm_pim; 

% auto volumes at ED and ES phases
compare_result.auto_esv_pic = lv_auto.esv_pic;
compare_result.auto_esv_pim = lv_auto.esv_pim;
compare_result.auto_edv_pic = lv_auto.edv_pic;
compare_result.auto_edv_pim = lv_auto.edv_pim;
compare_result.auto_sv_pic = lv_auto.sv_pic;
compare_result.auto_sv_pim = lv_auto.sv_pim;

% manual volume
compare_result.manual_esv_pic = lv_manual.esv_pic;
compare_result.manual_esv_pim = lv_manual.esv_pim;
compare_result.manual_edv_pic = lv_manual.edv_pic;
compare_result.manual_edv_pim = lv_manual.edv_pim;
compare_result.manual_sv_pic = lv_manual.sv_pic;
compare_result.manual_sv_pim = lv_manual.sv_pim;

compare_result.manual_ef_pic = lv_manual.ef_pic; %manual's ef
compare_result.manual_ef_pim = lv_manual.ef_pim; 
compare_result.manual_lvm_pic = lv_manual.lvm_pic; %manual's lv mass
compare_result.manual_lvm_pim = lv_manual.lvm_pim; 

compare_result.avg_dist_i = mean(avg_dist_i(auto_good_ix_i)); %average distance
compare_result.avg_dist_o = mean(avg_dist_o(auto_good_ix_o));
compare_result.avg_dm_i = mean(dm_i(auto_good_ix_i)); %average dice metric
compare_result.avg_dm_o = mean(dm_o(auto_good_ix_o));


end %compare contour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sindex = add_zero_index(index,digit_length)
    %add zero before index with total length of digit_length, for instance, change 20 to 0020.
    sindex=int2str(index);
    while length(sindex) < digit_length
           sindex = ['0', sindex]; %#ok<AGROW>
    end
end 

function full_contour_filename = get_contour_filename(contour_path, name_prefix, sindex, contour_mode, seg_mode)
    %gef contour filename with full path
    full_contour_filename = [contour_path filesep name_prefix, sindex, '-',contour_mode,'-',seg_mode, '.txt'];
end

function varargout = calc_dist(target_xy,reference_xy,contour_filename,para)
%calc perpendicular distance between target and reference contours

    %plot contours
    if para.save_figure
        %create distance figure folder
        if (para.auto_based_noraml == 1)
            str = strrep(contour_filename, para.auto_contour_foldername, [para.distance_figure_foldername '_auto_based_normal']);
        else
            str = strrep(contour_filename, para.auto_contour_foldername, [para.distance_figure_foldername '_manual_based_normal']);
        end
        
        [pathstr, name, ext] = fileparts(str);
        if (~exist(pathstr,'dir'))
            mkdir(pathstr)
        end
        
        %read dicom image
        str_dicom = strrep(contour_filename, [para.auto_contour_foldername filesep  para.auto_contour_subfoldername],para.dicom_foldername);
        dicom_filename = [str_dicom(1:end-18),'.dcm']; %IM-0001-XXXX-icontour-auto.txt
        raw_image = dicomread(dicom_filename);
        %show box
        min_showbox_x = round(min([target_xy(:,1); reference_xy(:,1)]));
        max_showbox_x = round(max([target_xy(:,1); reference_xy(:,1)]));
        min_showbox_y = round(min([target_xy(:,2); reference_xy(:,2)]));
        max_showbox_y = round(max([target_xy(:,2); reference_xy(:,2)]));
        %validate
        if min_showbox_x <1
           min_showbox_x = 1;
        end
        
        if min_showbox_y <1
           min_showbox_y = 1;
        end
                
        showbox_offset= 3;
        if ((min_showbox_x(1)-showbox_offset)<1 || (min_showbox_y(1)-showbox_offset)<1 || (max_showbox_x(1)+showbox_offset)>size(raw_image,2) || (max_showbox_y(1)+showbox_offset)>size(raw_image,1))
            showbox_offset = 0;
        end
        %crop image
        croped_image = raw_image(min_showbox_y(1)-showbox_offset:max_showbox_y(1)+showbox_offset,min_showbox_x(1)-showbox_offset:max_showbox_x(1)+showbox_offset);
        
        %show image
        h = figure;
        set(h,'Visible','off');
        axis image; clf;
        imshow(croped_image,[],'InitialMagnification',1000,'Border','tight')
        
        %plot contours
        
         if (para.auto_based_noraml == 1) %calc normal based on auto contours, manual constour is the reference contour
             hold on; plot(target_xy(:,1)-min_showbox_x(1)+showbox_offset,target_xy(:,2)-min_showbox_y(1)+showbox_offset,'r.-');
             hold on; plot(reference_xy(:,1)-min_showbox_x(1)+showbox_offset,reference_xy(:,2)-min_showbox_y(1)+showbox_offset,'g.-');
         else
             hold on; plot(reference_xy(:,1)-min_showbox_x(1)+showbox_offset,reference_xy(:,2)-min_showbox_y(1)+showbox_offset,'r.-');
             hold on; plot(target_xy(:,1)-min_showbox_x(1)+showbox_offset,target_xy(:,2)-min_showbox_y(1)+showbox_offset,'g.-');
         end
                
     end
    
    %check closed contours
    if target_xy(1,1) == target_xy(end,1) && target_xy(1,2) == target_xy(end,2)
        target_xy = target_xy(1:end-1,:);
    end
    if reference_xy(1,1) == reference_xy(end,1) && reference_xy(1,2) == reference_xy(end,2)
        reference_xy = reference_xy(1:end-1,:);
    end

    %for finding a subset of reference that correspond to the current target
    target_mask = poly2mask (target_xy(:,1),target_xy(:,2),double(para.width),double(para.height));
    reference_mask = poly2mask (reference_xy(:,1),reference_xy(:,2),double(para.width),double(para.height));
    
    %centroid of refernece mask
    reference_cen = regionprops(bwlabel(reference_mask),'Centroid');
    reference_cen_x = round(reference_cen.Centroid(1));
    reference_cen_y = round(reference_cen.Centroid(2));
    
    %check if target mask's centroid is in refernece's mask
    if  target_mask(reference_cen_y,reference_cen_x) ==1
        target_cen = reference_cen;
    else
        target_cen = regionprops(bwlabel(target_mask),'Centroid');
    end

    %calc degrees between contours points and centroid
    target_degree = atan2( target_xy(:,2)-target_cen.Centroid(2),target_xy(:,1)-target_cen.Centroid(1))*180/pi;
    reference_degree = atan2( reference_xy(:,2)-reference_cen.Centroid(2),reference_xy(:,1)-reference_cen.Centroid(1))*180/pi;
        
    %initialize
    distance=zeros(size(target_xy,1),1);
    
    for idx = 1: size(target_xy,1)
        %Get 3 points,[left current right]
        if idx == 1
            target_xy_3 = [target_xy(end,:); target_xy(idx:idx+1,:)];
        elseif (idx == size(target_xy,1))
            target_xy_3 = [target_xy(end-1:end,:); target_xy(1,:)];
        else
            target_xy_3 = target_xy(idx-1:idx+1,:);
        end

        %current target to reference distance
        target_curr_x = target_xy_3(2,1);
        target_curr_y = target_xy_3(2,2);

        %find a subset of reference that correspond to the current target
        if ( abs(target_degree(idx)) >= 150 )
            reference_xy_subset = reference_xy(find(abs(reference_degree)>=130),:);
        else
            reference_xy_subset = reference_xy(find(abs(reference_degree-target_degree(idx))<50),:);
        end
        
        %plot reference_xy_subset
        %hold on; plot(reference_xy_subset(:,1)-min_showbox_x(1)+showbox_offset, reference_xy_subset(:,2)-min_showbox_y(1)+showbox_offset ,'r*')
        
        %fit line
        target_X = target_xy_3(:,1);
        target_Y = target_xy_3(:,2);
        
        %3 points are the same
        if ( sum(abs(target_X - target_X(2))) + sum(abs(target_Y - target_Y(2))) == 0 )
            continue;
        end
        
        target_LINE = [target_X ones(size(target_X))] \ target_Y; %left division-least squares,  %y=LINE(1)*x+LINE(2)
        target_LINE_inv = [target_Y ones(size(target_Y))] \ target_X; %left division-least squares,  %x=LINE(1)*y+LINE(2)  Reverse the coordinates.
        
        %if the line is very close to parallel with the Y axis, and the residuals for the fit
        %improve by reversing the coordinates, then use the line found by
        %reversing the coordinates.
        if ( abs(target_LINE_inv(1))< 0.1 && sum(abs([target_Y ones(size(target_Y))] * target_LINE_inv - target_X)) < sum(abs([target_X ones(size(target_X))] * target_LINE - target_Y)))
            target_LINE = 1./(target_LINE_inv + 0.00001) ;
            %hold on; plot(target_X- min_showbox_x(1) + showbox_offset,target_Y - min_showbox_y(1) + showbox_offset,'Y*')
        end
        
        %Determine the line normal to the above tangent line (target_LINE).
        %Normal line: Ax+By+C=0;
        if abs(target_LINE(1)) < 0.01 %fit line parallel to X axis, normal parallel to Y axis
            A_N = 1;
            B_N = 0;
            C_N = -target_curr_x;
            reference_xy_subset_to_normal_distance = abs(reference_xy_subset(:,1)- target_curr_x); 
        elseif abs(target_LINE(1)) > 100 %fit line parallel to Y axis, normal parallel to X axis
            A_N = 0;
            B_N = 1;
            C_N = -target_curr_y;
            reference_xy_subset_to_normal_distance = abs(reference_xy_subset(:,2)- target_curr_y); 
        else
            A_N = -1/target_LINE(1);
            B_N = -1;
            C_N = target_curr_y - A_N * target_curr_x;
            reference_xy_subset_to_normal_distance = abs( (A_N*reference_xy_subset(:,1) + B_N*reference_xy_subset(:,2) + C_N))/sqrt(A_N^2 + B_N^2);%normal: y=N(1)*x+N(2);
        end

        %find reference point most close to normal
        [min_dist_subset, min_idx_subset] = min(reference_xy_subset_to_normal_distance);
      
        %if there are more than one min distance points, choose the one most close to current target
        min_idx_subset_1= find(reference_xy_subset_to_normal_distance <= 0.5); 
        
        if (length(min_idx_subset_1)>=2)
             min_dist_subset_to_target_curr_distance = sqrt((target_curr_x -  reference_xy_subset(min_idx_subset_1,1)).^2 + (target_curr_y - reference_xy_subset(min_idx_subset_1,2)).^2);
             [min_dist_subset, min_idx_subset_2] = min(min_dist_subset_to_target_curr_distance);
             reference_to_min_dist_subset_distance = sqrt((reference_xy(:,1) -  reference_xy_subset(min_idx_subset_1(min_idx_subset_2),1)).^2 + (reference_xy(:,2) - reference_xy_subset(min_idx_subset_1(min_idx_subset_2),2)).^2);
        else
             reference_to_min_dist_subset_distance = sqrt((reference_xy(:,1) -  reference_xy_subset(min_idx_subset,1)).^2 + (reference_xy(:,2) - reference_xy_subset(min_idx_subset,2)).^2);
        end
         
        min_idx_temp = find(reference_to_min_dist_subset_distance == 0);
        min_idx = min_idx_temp(1);

        %Get 3 points of reference,[left most-close-to-normal right]
        if min_idx == 1
            reference_xy_3 = [reference_xy(end,:); reference_xy(min_idx:min_idx+1,:)];
        elseif (min_idx == size(reference_xy,1))
            reference_xy_3 = [reference_xy(end-1:end,:); reference_xy(1,:)];
        else
            reference_xy_3 = reference_xy(min_idx-1:min_idx+1,:);
        end

        %fit reference line: Ax+By+C=0; 
        reference_X = reference_xy_3(:,1);
        reference_Y = reference_xy_3(:,2);
        reference_LINE = [reference_X ones(size(reference_X))] \ reference_Y; %left division-least squares,         %y=LINE(1)*x+LINE(2)
      
        reference_LINE_inv = [reference_Y ones(size(reference_Y))] \ reference_X; %left division-least squares,  %y=LINE(1)*x+LINE(2)
      
        %fit line parallel to Y axis,
        if ( abs(reference_LINE_inv(1))< 0.1 && sum(abs([reference_Y ones(size(reference_Y))] * reference_LINE_inv - reference_X)) < sum(abs([reference_X ones(size(reference_X))] * reference_LINE - reference_Y)))
            reference_LINE = 1./(reference_LINE_inv + 0.00001) ;
            %hold on; plot(target_X- min_showbox_x(1) + showbox_offset,target_Y - min_showbox_y(1) + showbox_offset,'Y*')
        end
        
        if abs(reference_LINE(1)) < 0.01 % reference line parallel to X axis
            A_R = 0;
            B_R = 1;
            C_R = -reference_Y(2);
        elseif abs(reference_LINE(1)) >100  % reference line parallel to Y axis
            A_R = 1;
            B_R = 0;
            C_R = -reference_X(2);
        else
            A_R = reference_LINE(1);
            B_R = -1;
            C_R = reference_Y(2) - A_R * reference_X(2);
        end
        
        %intersection of normal and reference line
        AA=[A_N B_N; A_R B_R];
        BB =[-C_N; -C_R];
        XX =AA\BB; %  AA is a square matrix, AA\BB is roughly the same as inv(AA)*BB

        %validate intersection point
        XX_reference_distance = sqrt((XX(1) -  reference_xy_subset(:,1)).^2 + (XX(2) - reference_xy_subset(:,2)).^2  );
        if min(XX_reference_distance)<1.5
            %distance
            distance(idx) = sqrt((target_curr_x - XX(1))^2 + (target_curr_y - XX(2))^2);
       
            if para.save_figure
                %plot intersection point
                hold on; plot(XX(1)-min_showbox_x(1)+showbox_offset,XX(2)-min_showbox_y(1)+showbox_offset,'m*') 
                %plot distance line
                hold on; plot([target_curr_x - min_showbox_x(1) + showbox_offset XX(1) - min_showbox_x(1) + showbox_offset], [target_curr_y - min_showbox_y(1) + showbox_offset XX(2) - min_showbox_y(1) + showbox_offset],'b'); 
            end
        else
            distance(idx) = -999;
        end
    end
    
    %in mm
    distance = distance * para.pixel_spacing(1);
    distance_effective = distance(distance>0);
    varargout{1} =mean(distance_effective); 
    
    if para.save_figure
        %plot minimum and maximum distance point
        min_dist_idx = find(distance == min(distance_effective));
        max_dist_idx = find(distance == max(distance_effective));
        hold on; plot(target_xy(min_dist_idx(1),1)-min_showbox_x(1)+showbox_offset,target_xy(min_dist_idx(1),2)-min_showbox_y(1)+showbox_offset,'yo')
        hold on; plot(target_xy(max_dist_idx(1),1)-min_showbox_x(1)+showbox_offset,target_xy(max_dist_idx(1),2)-min_showbox_y(1)+showbox_offset,'yd')
        
        %text minimum and maximum distance
        text(target_xy(min_dist_idx(1),1)-min_showbox_x(1)+showbox_offset+0.5,target_xy(min_dist_idx(1),2)-min_showbox_y(1)+showbox_offset+0.5,[num2str(distance(min_dist_idx(1))),'(mm)'],'color','g','FontSize',11) %+0.5, not polt overlap
        text(target_xy(max_dist_idx(1),1)-min_showbox_x(1)+showbox_offset+0.5,target_xy(max_dist_idx(1),2)-min_showbox_y(1)+showbox_offset+0.5,[num2str(distance(max_dist_idx(1))),'(mm)'],'color','g','FontSize',11)
        
        %text mean and std of distances
        text(2,2,['Mean:' num2str(mean(distance_effective)) '(mm)' ', Std:' num2str(std(distance_effective)) '(mm)'],'color','g','FontSize',11);
        
        %legend
        legend('target','reference','intersection','distance line')
        
        %save figure as png
        saveas(h,  [str(1:end-9), '.png']); %IM-0001-XXXX-icontour-auto.txt
        close(h)
    end
end

function dm = calc_dm(autoPoints,manualPoints,para)
    %calc dice metric
    auto_mask = poly2mask (autoPoints(:,1),autoPoints(:,2),double(para.width),double(para.height));
    manual_mask = poly2mask (manualPoints(:,1),manualPoints(:,2),double(para.width),double(para.height));

    auto_size = sum(auto_mask(:)>0);
    manual_size = sum(manual_mask(:)>0);
    intersect_size = sum((auto_mask(:) + manual_mask(:))==2);
    dm = 2 * intersect_size / (auto_size + manual_size);
end

function lv = calc_clinical_para(dicomPath, contour_path, para,varargin)
    %calc ejection fraction and lv mass
       
    [dicomCount, sliceCount] = dicom_counter(dicomPath,para.phase_number);
    contTable_i = zeros(sliceCount, para.phase_number);
    contTable_o = zeros(sliceCount, para.phase_number);
    contTable_p1 = zeros(sliceCount, para.phase_number);
    contTable_p2 = zeros(sliceCount, para.phase_number);
    
    if length(varargin) > 3
        startingSlice = varargin{1};
        endingSlice = varargin{2};
        systolePhase = varargin{3};
        diastolePhase = varargin{4};
    elseif length(varargin) > 1
        startingSlice = varargin{1};
        endingSlice = varargin{2};
        systolePhase = 0;
        diastolePhase = 0;
    else 
        startingSlice = 1;
        endingSlice = sliceCount;
        systolePhase = 0;
        diastolePhase = 0;    
    end    
    
    if startingSlice < 1 || startingSlice > sliceCount
        startingSlice = 1;
    end

    if endingSlice < 1 || endingSlice > sliceCount
        endingSlice = sliceCount;
    end
    
    if startingSlice > endingSlice
        endTemp = startingSlice;
        startingSlice = endingSlice;
        endingSlice = endTemp;
    end   
    
    %list all contours
    contour_list = dir([contour_path filesep '*contour*.txt']);
    
    for contourIdx = 1:length(contour_list)
        contourFile = contour_list(contourIdx).name;

        %inside contour
        if (strcmp(contourFile(14),'i'))
            
            %manual contours
            if para.auto == false
               good_icontour = true;
            end
            %auto contours
            if para.auto == true
               good_icontour = any(ismember(para.auto_good_index_i, str2num(contourFile(9:12)))); 
            end
        
            if ~good_icontour
                continue;
            end
            
            match = regexp((contourFile),'-','start');
            imNum = contourFile(match(2) + 1:match(3) - 1);
            currSlice = getSlice(str2double(imNum),para.phase_number);
            currPhase = getPhase(str2double(imNum),para.phase_number);
            full_contour_filename = [contour_path filesep contourFile];
           try
               xy = load(full_contour_filename);
               area = polyarea(xy(:,1),xy(:,2));
               area_mm = area* para.pixel_spacing(1)^2;
               vol_cm3 = area_mm * (para.thickness + para.gap) *(1/10)^3; %change to cm3
               contTable_i(currSlice, currPhase) = vol_cm3; 
           catch
               s = lasterror;
               disp(s.message);
               %return;
            end
            
        end
        
        %outside contour
        if (strcmp(contourFile(14),'o') )
            if para.auto == false
                good_ocontour = true;
            end
            if para.auto == true
               good_ocontour = any(ismember(para.auto_good_index_o, str2num(contourFile(9:12)))); 
            end
            if ~good_ocontour
                continue;
            end
            
            match = regexp((contourFile),'-','start');
            imNum = contourFile(match(2) + 1:match(3) - 1);
            currSlice = getSlice(str2double(imNum),para.phase_number);
            currPhase = getPhase(str2double(imNum),para.phase_number);
            full_contour_filename = [contour_path filesep contourFile];
            try
               xy = load(full_contour_filename);
               area = polyarea(xy(:,1),xy(:,2));
               area_mm = area* para.pixel_spacing(1)^2;
               vol_cm3 = area_mm * (para.thickness + para.gap) *(1/10)^3;%change to cm3
               contTable_o(currSlice, currPhase) = vol_cm3; 
           catch
               s = lasterror;
               disp(s.message);
               %return;
            end
            
        end
        
        %pap contour
        if (strcmpi(contourFile(14), 'p') )
            if para.auto == false
                good_pcontour = true;
            end
            if para.auto == true
                good_pcontour = any(ismember([ para.auto_good_index_i; para.auto_good_index_o] , str2num(contourFile(9:12))));
            end
            if ~good_pcontour
                continue;
             end
            
            match = regexp((contourFile),'-','start');
            imNum = contourFile(match(2) + 1:match(3) - 1);
            currSlice = getSlice(str2double(imNum),para.phase_number);
            currPhase = getPhase(str2double(imNum),para.phase_number);
            full_contour_filename = [contour_path filesep contourFile];
            try
               xy = load(full_contour_filename);
               area = polyarea(xy(:,1),xy(:,2));
               area_mm = area* para.pixel_spacing(1)^2;
               vol_cm3 = area_mm * (para.thickness + para.gap) *(1/10)^3;%change to cm3
               
               if ( strcmpi(contourFile(15), '1') )
                    contTable_p1(currSlice, currPhase) = vol_cm3; 
               end
               if ( strcmpi(contourFile(15), '2') )
                    contTable_p2(currSlice, currPhase) = vol_cm3; 
               end
            catch
               s = lasterror;
               disp(s.message);
               %return;
            end
        end %end if
        
    end     
  
    maxPhase = length(contTable_i(1,:));
    if systolePhase < 1 || systolePhase > maxPhase || diastolePhase < 1 || diastolePhase > maxPhase
        [systolePhase, diastolePhase] = esedDetermine(contTable_i);
    end
    
    [contTableChecked_i, zeroed, slicesExcluded, slicesIncluded] = efCheck(contTable_i, systolePhase, diastolePhase, startingSlice, endingSlice);
    
    %calc ef - use only inside contour info.
    esvol_pic = sum(contTableChecked_i(:, systolePhase)); %'_pic' means pap included in the LV cavity
    esvol_pim = sum(contTableChecked_i(:, systolePhase))- sum(contTable_p1(slicesIncluded,systolePhase))- sum(contTable_p2(slicesIncluded,systolePhase));%'_pim' means pap included in the myocardium
    edvol_pic = sum(contTableChecked_i(:, diastolePhase));
    edvol_pim = sum(contTableChecked_i(:, diastolePhase)) - sum(contTable_p1(slicesIncluded,diastolePhase))- sum(contTable_p2(slicesIncluded,diastolePhase));
    strokeVol_pic = edvol_pic - esvol_pic;
    strokeVol_pim = edvol_pim - esvol_pim;
    ef_pic = strokeVol_pic / edvol_pic * 100;
    ef_pim = strokeVol_pim / edvol_pim * 100;
    
    lv.esv_pic = esvol_pic;
    lv.esv_pim = esvol_pim;
    lv.edv_pic = edvol_pic;
    lv.edv_pim = edvol_pim;
    lv.sv_pic = strokeVol_pic;
    lv.sv_pim = strokeVol_pim;
    lv.ef_pic = ef_pic;
    lv.ef_pim = ef_pim;
    lv.es = systolePhase;
    lv.ed = diastolePhase;
    lv.zeroed = zeroed;
    lv.slice_excluded = slicesExcluded;
   
    %calc lv mass 
    combined_table(:,1) = contTable_i(:,diastolePhase);
    combined_table(:,2) = contTable_o(:,diastolePhase);
    common_slice = ((combined_table(:,1) ~= 0) + (combined_table(:,2) ~= 0)) == 2 ;
    
    edv_i_pic = sum(combined_table(common_slice, 1));
    edv_o_pic = sum(combined_table(common_slice, 2));
    
    lvm_pic = (edv_o_pic - edv_i_pic) * 1.05; %1.05 (g/cm3)
    lvm_pim = (edv_o_pic - edv_i_pic + sum(contTable_p1(common_slice,diastolePhase)) + sum(contTable_p2(common_slice,diastolePhase)) ) * 1.05; 
    
    lv.lvm_pic = lvm_pic;
    lv.lvm_pim = lvm_pim;
end

function [systolePhase, diastolePhase] = esedDetermine(contTable)
    %determin ES and ED phase 

    contVolSum = sum(contTable);
    minVol = 9999; 
    maxVol = 0;
    
    for x = 1:length(contVolSum)
        valueUnderExamination = contVolSum(x);
        
        if valueUnderExamination > 0
            if minVol > valueUnderExamination
                minVol = contVolSum(x);
                systolePhase = x;
            end

            if maxVol < valueUnderExamination
                maxVol = contVolSum(x);
                diastolePhase = x;
            end
        end
    end
end

function sliceNum = getSlice(imageNum,phase_number)
    % calc slice number from the image number
     sliceNum = ceil(imageNum/phase_number);
end

function phaseNum = getPhase(imageNum,phase_number)
    %calc phase number from image number
    remainder = mod(imageNum,phase_number);
    if (remainder ~= 0)
        phaseNum = remainder;
    else
        phaseNum = phase_number;
    end 
end

function [volTableCropped, zeroed, slicesExcluded,slicesIncluded] = efCheck(volTable, systolePhase, diastolePhase, startSlice, endSlice)
    % modifies the contour table that is passed in by...
    % 1. removing the slices that are to be excluded using startSlice and
    % endSlice by zeroing all those values
    % 2. removing ED/ES pairs that are missing the other value.  
    % that is, if there is a zero in the systolic phase but not in the
    % diastolic phase in that slice, zero the other
    
    % 1. removing the non-included slices
    volTableCropped = zeros(length(volTable(:,1)), length(volTable(1,:)));
    volTableCropped(startSlice:endSlice, :) = volTable(startSlice:endSlice, :);
    
    % 2. finding the ES/ED slices with missing values
    sysZeros = find(volTable(:, systolePhase) == 0);
    diaZeros = find(volTable(:, diastolePhase) == 0);
    
    if ~isempty(sysZeros)
        volTableCropped(sysZeros, :) = 0;
    end
    
    if ~isempty(diaZeros)
        volTableCropped(diaZeros, :) = 0;
    end
    
    slicesExcluded = find(volTableCropped(:, systolePhase) == 0);
    slicesIncluded = find(volTableCropped(:, systolePhase) ~= 0);
    
    zeroed = length(slicesExcluded);
end

function [dicom_count, slice_count] = dicom_counter(dicom_path, phase_number)
    %calc dicom image number and slice number
    dicom_count = 0;
    slice_count = 0;
    
    if ~exist(dicom_path,'dir')
        %return
    end
    
    dicom_files = dir([dicom_path filesep '*.dcm']);
    if isempty(dicom_files)
        disp([dicom_path ' :NO DICOM files!'])
        %return;
    else 
        dicom_count = length(dicom_files);
        slice_count = ceil(dicom_count / phase_number);
    end
end
