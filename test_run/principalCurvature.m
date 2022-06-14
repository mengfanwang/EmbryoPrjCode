clear;
im = im2double(tiffRead('MAX_fish13_after_treatment_0923.tif',16));

im = imgaussfilt(im,1.25);
imagesc(im);
[lx, ly] = gradient(im);
[lxx,lyx] = gradient(lx);
[lxy, lyy] = gradient(ly);

eig1 = zeros(size(im));
eig2 = zeros(size(im));
u1 = zeros(size(im));
v1 = zeros(size(im));
u2 = zeros(size(im));
v2 = zeros(size(im));
for i  =1: length(im(:))
    MM = [lxx(i) lxy(i);lyx(i) lyy(i)];
    [Evec,Eval] = eig(MM);
    dEval = diag(Eval);
    [c , id] = sort(dEval);
    eig1(i) = dEval(id(1));
    eig2(i) = dEval(id(2));
    u1(i) = Evec(1,id(1));
    v1(i) = Evec(2,id(1));
    u2(i) = Evec(1,id(2));
    v2(i) = Evec(2,id(2));
end

[x,y] = meshgrid(1:1767,1:399);
figure(1)
imagesc(eig1);colorbar
hold on
quiver(x, y,u1,v1)



[x,y] = meshgrid(1:1767,1:399);
figure(2)
imagesc(eig2);colorbar
hold on
quiver(x, y,u2,v2)



imagesc(im)


%% %%%%%%%%%%%%%%%3D%%%%%%%%%%%%%%%%%%%%%%%%%%
im = scale_image(crop_embryo_vid(:,:,:,1),0,1);
% im = im(:,:,67:123);
sigma = 2;
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(im,sigma);

if(sigma>0)
    % Correct for scaling
    c=(sigma^2);
    Dxx = c*Dxx; Dxy = c*Dxy;
    Dxz = c*Dxz; Dyy = c*Dyy;
    Dyz = c*Dyz; Dzz = c*Dzz;
end
eig1 = zeros(size(im));
eig2 = zeros(size(im));
eig3 = zeros(size(im));

impr = zeros(size(im));
u1 = zeros(size(im));
v1 = zeros(size(im));
w1 = zeros(size(im));
u2 = zeros(size(im));
v2 = zeros(size(im));
w2 = zeros(size(im));

u3 = zeros(size(im));
v3 = zeros(size(im));
w3 = zeros(size(im));


parfor i  =1: length(impr(:))
    MM = [Dxx(i), Dxy(i), Dxz(i);Dxy(i), Dyy(i), Dyz(i);Dxz(i), Dyz(i), Dzz(i)];
    
    [Evec,Eval] = eig(MM);
    %c = sort(diag(Eval));
    dEval = diag(Eval);
    [c , id] = sort(dEval,'descend');
    eig1(i) = dEval(id(1));
    eig2(i) = dEval(id(2));
    eig3(i) = dEval(id(3));
    u1(i) = Evec(1,id(1));
    v1(i) = Evec(2,id(1));
    w1(i) = Evec(3,id(1));
    u2(i) = Evec(1,id(2));
    v2(i) = Evec(2,id(2));
    w2(i) = Evec(3,id(2));
    u3(i) = Evec(1,id(3));
    v3(i) = Evec(2,id(3));
    w3(i) = Evec(3,id(3));
    
end
[lx, ly, lz] = size(im);
figure(4)
imagesc(im(:,:,17))
colorbar
[x,y] = meshgrid(1:ly,1:lx);
u11 = u1(:,:,17);
v11 = v1(:,:,17);
w11 = w1(:,:,17);     
figure(1)
imagesc(eig1(:,:,17));colorbar
hold on
quiver(x, y,v11,u11)
title('\lambda_{1}')

[x,y] = meshgrid(1:ly,1:lx);
u22 = u2(:,:,17);
v22 = v2(:,:,17);
w22 = w2(:,:,17);
figure(6)
imagesc(eig2(:,:,17));colorbar
hold on
quiver(x,y,v22,u22);
title('\lambda_{2}')

[x,y] = meshgrid(1:ly,1:lx);
u22 = u2(:,:,17);
v22 = v2(:,:,17);
w22 = w2(:,:,17);
figure(2)
imagesc(im(:,:,17));colorbar
hold on
quiver(x,y,v22,u22,'r');
title('\lambda_{2} on original image')

[x,y] = meshgrid(1:ly,1:lx);
u33 = u3(:,:,17);
v33 = v3(:,:,17);
w33 = w3(:,:,17);
figure(3)
imagesc(eig3(:,:,17));colorbar
hold on
quiver(x,y,v33,u33);
title('\lambda_{3}')


% [x,y] = meshgrid(1:1036,1:565);
% u33 = u3(:,:,4);
% v33 = v3(:,:,4);
% w33 = w3(:,:,4);
% figure(7)
% imagesc(im(:,:,4));colorbar
% hold on
% quiver(x,y,v33,u33, 'r');
% title('\lambda_{3} on original image')


%%%%%%%%%%show the eig values on the yz dimension%%%%%%%%%%%%%%%%%%

figure(4)
imagesc(squeeze(im(:,377,:)))
colorbar


[x,y] = meshgrid(1:35,1:277);
u33 = u3(:,377,:);
u33 = squeeze(u33);
v33 = v3(:,377,:);
v33 = squeeze(v33);
w33 = w3(:,377,:);
w33 = squeeze(w33);
figure(3)
imagesc(squeeze(eig3(:,377,:)));colorbar
hold on
quiver(x,y,v33,u33);
title('\lambda_{3}')


[x,y] = meshgrid(1:35,1:277);
u22 = u2(:,377,:);
u22 = squeeze(u22);
v22 = v2(:,377,:);
v22 = squeeze(v22);
w22 = w2(:,377,:);
w22 = squeeze(w22);
figure(2)
imagesc(squeeze(eig2(:,377,:)));colorbar
hold on
quiver(x,y,v22,u22);
title('\lambda_{2}')


[x,y] = meshgrid(1:35,1:277);
u11 = u1(:,377,:);
u11 = squeeze(u11);
v11 = v1(:,377,:);
v11 = squeeze(v11);
w11 = w3(:,377,:);
w11 = squeeze(w11);
figure(1)
imagesc(squeeze(eig1(:,377,:)));colorbar
hold on
quiver(x,y,v11,u11);
title('\lambda_{1}')









% u11 = u1(1:10:138,1:10:884,1:10:125);
% v11 = v1(1:10:138,1:10:884,1:10:125);
% w11 = w1(1:10:138,1:10:884,1:10:125);
impr = zeros(size(im));
alpha = 0.25;
for j = 1: length(impr(:))
    if(eig3(j)< 0 && eig2(j) <0)
        impr(j) = abs(eig2(j)) + abs(eig3(j));
    %elseif(eig2(j) <0 && eig3(j)>0 && eig3(j)<abs(eig2(j))/alpha)
        %impr(j) = abs(eig2(j)) - alpha* eig3(j);
    elseif(eig3(j) <0 && eig2(j) >0)
        impr(j) = abs(eig3(j)) - alpha*eig2(j);
    else
        impr(j) = 0;
    end
end
% impr = sqrt(eig1x.*eig2x);
impr = impr./max(impr(:));
impr2 = abs(impr);
impr2(impr2(:) < 0.1) = 0;
impr2(impr2(:)~=0) = 1;
volumeViewer(impr2)

%%%%%%%%%%%%%%%If we select the convex regions only

impr = zeros(size(im));
for j = 1: length(impr(:))
    if(eig3(j)<0 && eig2(j)<0)
        impr(j) = abs(eig2(j)) + abs(eig3(j));
    end
end

impr = impr./max(impr(:));
impr2 = abs(impr);
impr2(impr2(:) < 0.1) = 0;
impr2(impr2(:)~=0) = 1;
volumeViewer(impr2)

% 
% impr = eig2;
% impr = impr./max(impr(:));
% impr2 = abs(impr);
% impr2(impr2(:) < 0.05) = 0;


tifwrite(im2double(impr2), 'impr2_fish13_1016_large_finerSegmentation');


