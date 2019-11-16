%% DualTree Complexe Wavelet 

%% Reading the gif image to RGB
[A_CT,map_CT]= imread('023.gif',1);
RGB_CT = ind2rgb(A_CT,map_CT);
[A_S,map_S]   = imread('023-1.gif',1);
RGB_S = ind2rgb(A_S,map_S);

%% Forward Transform of Images
% CT-images
impl_CT  = dddtree2('cplxdt', im2unit8(RGB_CT(:,:,1)) ,5,'dtf3');
impl1_CT = dddtree2('cplxdt', im2unit8(RGB_CT(:,:,2)) ,5,'dtf3');
impl2_CT = dddtree2('cplxdt',  im2unit8(RGB_CT(:,:,3)) ,5,'dtf3');
% SPET Images
impl_S  = dddtree2('cplxdt', im2unit8(RGB_S(:,:,1)) ,5,'dtf3');
impl1_S = dddtree2('cplxdt', im2unit8(RGB_S(:,:,2)) ,5,'dtf3');
impl2_S = dddtree2('cplxdt', im2unit8(RGB_S(:,:,3)) ,5,'dtf3');

%  % Highpass Coefficients
% for j = 1:5
%     for d = 1:3
%         for k = 1:2
%             for m = 1:2
% F_Coe = (impl_CT.cfs{j}(:,:,d,k,m)+impl_S.cfs{j}(:,:,d,k,m))/2;
% F_Coe1 = (impl1_CT.cfs{j}(:,:,d,k,m)+impl1_S.cfs{j}(:,:,d,k,m))/2;
% F_Coe2 = (impl2_CT.cfs{j}(:,:,d,k,m)+impl2_S.cfs{j}(:,:,d,k,m))/2;
% impl_CT.cfs{j}(:,:,d,k,m)=F_Coe;
% impl1_CT.cfs{j}(:,:,d,k,m)=F_Coe1;
% impl2_CT.cfs{j}(:,:,d,k,m)=F_Coe2;
%             end
%         end
%     end
% end
%% Lowpass Coefficients 
a = 1
b = 1
for j = 1:5
    for n = 1:8
        for g = 1:8
            for k = 1:2
                for m = 1:2
F_Coe = (a*impl_CT.cfs{j}(n,g,k,m)+b*impl_S.cfs{j}(n,g,k,m))/2;
F_Coe1 = (a*impl1_CT.cfs{j}(n,g,k,m)+b*impl1_S.cfs{j}(n,g,k,m))/2;
F_Coe2 = (a*impl2_CT.cfs{j}(n,g,k,m)+b*impl2_S.cfs{j}(n,g,k,m))/2;
impl_CT.cfs{j}(n,g,k,m)=F_Coe;
impl1_CT.cfs{j}(n,g,k,m)=F_Coe1;
impl2_CT.cfs{j}(n,g,k,m)=F_Coe2;
                end
            end 
        end 
    end 
end 
       
%% Inverse Transform
 dtImage1 = idddtree2(impl_CT);
 dtImage2 = idddtree2(impl1_CT);
 dtImage3 = idddtree2(impl2_CT);
%% regrouping of image channels
image(:,:,1) = dtImage1;
image(:,:,2) = dtImage2;
image(:,:,3) = dtImage3;
%% show
figure
imshow(image)

%% measure
e = entropy(image)
e1 = entropy(RGB_CT)
e2 = entropy(RGB_S)

