

function para=get_dicominfo(dicom_path)

% find dcm images
dicom_path1=[dicom_path,'/*.dcm'];
dicom_files=dir(dicom_path1);


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

end
