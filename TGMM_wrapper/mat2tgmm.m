clc;clear;close all;

load('E:\Embryo\TM0-49\track_v1\movieInfo.mat');
load('E:\Embryo\TM0-49\track_v1\track_refine_res.mat');
xml_folder = 'E:\Embryo\TM0-49\track_v1\tgmm_xml\';
% registration_flag = 0;

t = length(movieInfo.n_perframe);
start_ind = [0; cumsum(movieInfo.n_perframe)];
n_total = length(movieInfo.xCoord);
lineage_list = nan(n_total,1);
x_ind = movieInfo.xCoord;
y_ind = movieInfo.yCoord;
z_ind = movieInfo.zCoord;
% if registration_flag
%     % load trans_mat
%     load('E:\Embryo\registration_temporal_data\ImageReg.mat');
%     trans_before = [x_ind y_ind z_ind ones(length(x_ind),1)];
%     trans_after = zeros(length(x_ind),3);
%     for ii = 1:length(x_ind)
%         tt = find(ii > start_ind, 1, 'last');
%         trans_temp = trans_before(ii,:)*trans_mat{tt};
%         trans_after(ii,:) = trans_temp(1:3);
%     end
%     x_ind = trans_after(:,1) - bound.x_min;
%     y_ind = trans_after(:,2) - bound.y_min;
%     z_ind = trans_after(:,3) - bound.z_min;
% end
for ii = 1:length(movieInfo.tracks)
    lineage_list(movieInfo.tracks{ii}) = ii - 1;
end
for tt = 1:t
    tt
    docNode = com.mathworks.xml.XMLUtils.createDocument('document');
    root = docNode.getDocumentElement;
    for ii = start_ind(tt)+1:start_ind(tt+1) 
%     for ii = 1:1
        lineage = lineage_list(ii);
        if isnan(lineage)
            continue;
        end
        id = ii - (start_ind(tt)+1);
        parent = -1;
        if ~isempty(movieInfo.parents{ii})
            parent = movieInfo.parents{ii} - (start_ind(tt-1)+1);
        end
        
        gmm = docNode.createElement('GaussianMixtureModel');
        gmm.setAttribute('id',num2str(id));
        gmm.setAttribute('lineage',num2str(lineage));
        gmm.setAttribute('parent',num2str(parent));
        gmm.setAttribute('splitScore','3');
        gmm.setAttribute('scale','1 1 1');
        gmm.setAttribute('nu','100');
        gmm.setAttribute('beta','100');
        gmm.setAttribute('alpha','100');
        gmm.setAttribute('m',[num2str(x_ind(ii)) ' ' num2str(y_ind(ii)) ' ' num2str(z_ind(ii))]);
        gmm.setAttribute('W','0.01 0 0 0 0.01 0 0 0 0.01');
        gmm.setAttribute('svIdx',num2str(id));
        gmm.appendChild(docNode.createComment('.'));
        root.appendChild(gmm);
    end
    ind = num2str(10000+tt-1);
    ind = ind(2:end);
    xmlwrite([xml_folder 'GMEMfinalResult_frame' ind '.xml'], docNode);
end