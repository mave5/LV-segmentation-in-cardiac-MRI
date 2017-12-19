% save contours into txt files
function save_contours2(auto_mask,cotourfilename)
    
     % get manual contour file name
     %origStr=char(cont_name);

     % replace manual with auto in the name
     %modifiedStr = strrep(origStr, 'manual', 'auto');

     % convert mask to contours
     temp=(contourc(auto_mask,[0 0]))';
     auto_contour=temp(2:end,:);

     % save contours into txt files            
     filename=cotourfilename;
     save (filename,'auto_contour','-ASCII');

display('contour save into text file')
