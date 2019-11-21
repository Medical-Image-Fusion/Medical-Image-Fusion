% Get the list of images andcreate the folder for the output images
cd Images
if ~exist('Fused' , 'dir')
    mkdir('Fused')
end
cd MRI
MRIFiles = dir('*.gif');
cd ..
cd PET
PETFiles = dir('*.gif');
cd ..
cd ..

ifpm = zeros(length(MRIFiles), 9); %Create the basic ifmp array
fileOrder = repmat(' ', length(MRIFiles), 1); %Create the string array to know the order for which the imgs were loaded

parfor k = 1:length(MRIFiles)
  baseFileName = MRIFiles(k).name; %Get the name of one of the images on the list
  fileOrder(k) = MRIFiles(k).name; %Save the name on the array
  fullFileNameMRI = fullfile(MRIFiles(k).folder, baseFileName); %Get the file path for the MRI image
  fprintf(1, 'Now reading %s\n', fullFileNameMRI);
  fullFileNamePET = fullfile(PETFiles(k).folder, baseFileName); %Get the file path for the PET image
  fprintf(1, 'Now reading %s\n', fullFileNamePET);
  [A,map] = imread(fullFileNameMRI,1); %Load the RMI img
  A_RGB = ind2rgb(A,map);
  [B,map] = imread(fullFileNamePET,1); %Load the PET img
  B_RGB = ind2rgb(B,map);
  [ifpm_k, fused_img] = MWT_PSO(A_RGB, B_RGB); %Run the PSO function and return the fused img + the coeff used
  ifpm(k,:) = ifpm_k; %Assign the coeff to the array
  
  name = baseFileName(1:end-4); % Get the name of the img
  imwrite(fused_img, strcat('Images/Fused/Fused_', name,'.png')) %Save the fused img
  parsave(ifpm_k, name); % Save the .mat with the coefficients and ifpm

%   imageArray = imread(fullFileName);
%   imshow(imageArray);  % Display image.
%   drawnow; % Force display to update immediately.
end
save('image_fusion') % Save the whole workspace, just in case

function parsave(ifpm_k, name)
    save(strcat('Images/Fused/Fused_', name,'.mat'), 'ifpm_k');
end