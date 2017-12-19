function varargout = calc_dist(target_xy,reference_xy,para)
%calc perpendicular distance between target and reference contours

    
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
    reference_cen = regionprops(logical(reference_mask),'Centroid');
    reference_cen_x = round(reference_cen.Centroid(1));
    reference_cen_y = round(reference_cen.Centroid(2));
    
    %check if target mask's centroid is in refernece's mask
    if  target_mask(reference_cen_y,reference_cen_x) ==1
        target_cen = reference_cen;
    else
        %target_mask=clean_segs(target_mask);
        target_cen = regionprops(logical(target_mask),'Centroid');
        
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
       
        else
            distance(idx) = -999;
        end
    end
    
    %in mm
    distance = distance * para.pixel_spacing(1);
    distance_effective = distance(distance>0);
    varargout{1} =mean(distance_effective); 
    
end
