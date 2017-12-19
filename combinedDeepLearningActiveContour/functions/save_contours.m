% save contours into txt files
function output=save_contours(masks,t_cont_names,slice_per_patient)

% get number of studies
num_studies=length(slice_per_patient);


for k=1:num_studies
    
    % get number of slices per patient
    spp=slice_per_patient(k);

    % make a new directory
    dirn=['matFiles/auto_contours/auto',num2str(k)];
    mkdir(dirn);
    
    for k2=1:spp
        
         % get manual contour file name
         ind=k2+sum(slice_per_patient(1:k-1));
         origStr=t_cont_names(ind).name;

         % replace manual with auto in the name
         modifiedStr = strrep(origStr, 'manual', 'auto');

         % convert mask to contours
         temp=(contourc(masks(:,:,ind),[0 0]))';
         auto_contour=temp(2:end,:);
         
         % save contours into txt files         
         filename=[dirn,'/',modifiedStr];
         save (filename,'auto_contour','-ASCII');
    end  
end

output='done'
