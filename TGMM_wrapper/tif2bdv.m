clc;clear;close all;
% tic;hdf5info('E:\Embryo\TM0-49\H2BGFP_TM0-49.h5');toc
% read tif file
folder_path = 'E:\Embryo\registration_temporal_data\PointReg';
% folder_list = dir(folder_path);
% folder_list = folder_list(3:end);
data = zeros(243,380,79,40);
for tt = 0:39
%     file = dir([folder_path folder_list(tt).name]);
%     file = file(end);
    ind = num2str(100+tt);
    ind = ind(2:end);
    data(:,:,:,tt+1) = tifread([folder_path '\' ind '.tif']);
end

%%% parameter %%%
data_name = 'PointReg';
view_num = 1;
time_num = 40;
res_mat = [1 1 1; 2 2 2;];
res_mat = res_mat';
sub_mat = 16*ones(3,2);

% data transform
[x,y,z,t] = size(data);
data_t = zeros(y,x,z,t);
for tt = 1:t
    for zz = 1:z
        data_t(:,:,zz,tt) = data(:,:,zz,tt)';
    end
end
data = data_t;

% write xml file
docNode = com.mathworks.xml.XMLUtils.createDocument('SpimData');
root = docNode.getDocumentElement;
root.setAttribute('version','0.2');
basePath = docNode.createElement('BasePath');
basePath.setAttribute('type','relative');
basePath.appendChild(docNode.createTextNode('.'));
root.appendChild(basePath);

% sequence description
seqDes = docNode.createElement('SequenceDescription');

imLoader = docNode.createElement('ImageLoader');
imLoader.setAttribute('format','bdv.hdf5');
hdf5 = docNode.createElement('hdf5');
hdf5.setAttribute('type','relative');
hdf5.appendChild(docNode.createTextNode([data_name '.h5']));
imLoader.appendChild(hdf5);
seqDes.appendChild(imLoader);

viewSets = docNode.createElement('ViewSetups');
for ii = 0:view_num-1
    viewSet = docNode.createElement('ViewSetup');
    id = docNode.createElement('id');
    id.appendChild(docNode.createTextNode(num2str(ii)));   % id name
    viewSet.appendChild(id);
    viewSets.appendChild(viewSet);
end
seqDes.appendChild(viewSets);

timePoints = docNode.createElement('Timepoints');
timePoints.setAttribute('type','range');
first = docNode.createElement('first');
first.appendChild(docNode.createTextNode('0'));
timePoints.appendChild(first);
last = docNode.createElement('last');
last.appendChild(docNode.createTextNode(num2str(time_num-1)));
timePoints.appendChild(last);
seqDes.appendChild(timePoints);

root.appendChild(seqDes);

% view registration
viewRegs = docNode.createElement('ViewRegistrations');
for tt = 0:time_num-1
    for ii = 0:view_num-1 
        viewReg = docNode.createElement('ViewRegistration');
        viewReg.setAttribute('timepoint',num2str(tt));
        viewReg.setAttribute('setup',num2str(ii));
        viewTrans = docNode.createElement('ViewTransform');
        viewTrans.setAttribute('type','affine');
        affine = docNode.createElement('affine');
        affine.appendChild(docNode.createTextNode(...
            '1 0 0 0 0 1 0 0 0 0 1 0'));
        viewTrans.appendChild(affine);
        viewReg.appendChild(viewTrans);
        viewRegs.appendChild(viewReg);
    end
end
root.appendChild(viewRegs);
xmlwrite(['.\' data_name '.xml'],docNode);

% write hdf5 data
h5create([data_name '.h5'],'/s00/resolutions',size(res_mat));
h5write([data_name '.h5'], '/s00/resolutions', res_mat);
h5create([data_name '.h5'],'/s00/subdivisions',size(sub_mat));
h5write([data_name '.h5'], '/s00/subdivisions', sub_mat);
l_num = size(res_mat,2);
for tt = 0:time_num-1
    t_ind = num2str(100000+tt);
    t_ind = t_ind(2:end);
    s_ind = '00';
    for ll = 0:l_num-1
        l_ind = num2str(ll);
        data_temp = imresize3(data(:,:,:,tt+1),round(size(data(:,:,:,tt+1))./res_mat(:,ll+1)'));
        h5create([data_name '.h5'],['/t' t_ind '/s' s_ind '/' l_ind '/cells'],...
            size(data_temp),'Datatype','uint16','ChunkSize',sub_mat(:,ll+1)');
        h5write([data_name '.h5'],['/t' t_ind '/s' s_ind '/' l_ind '/cells'],data_temp);
    end
end
