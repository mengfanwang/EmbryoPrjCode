clc;clear;close all;

if isunix
    load('/work/Mengfan/Embryo/TM0-49/track_0.25/movieInfo_temp.mat');
    load('/work/Mengfan/Embryo/TM0-49/track_0.25/track_refine_res.mat');
    xml_folder = '/work/Mengfan/Embryo/TM0-49/tgmm_xml';
else
    load('E:\Embryo\TM0-49\track_v1\movieInfo.mat');
    load('E:\Embryo\TM0-49\track_v1\track_refine_res.mat');
    xml_folder = 'E:\Embryo\TM0-49\track_v1\tgmm_xml\';
end


t = length(movieInfo.n_perframe);
start_ind = [0; cumsum(movieInfo.n_perframe)];
n_total = length(movieInfo.xCoord);
lineage_list = nan(n_total,1);
x_ind = (movieInfo.xCoord - 1) * 2;
y_ind = (movieInfo.yCoord - 1) * 2;
z_ind = (movieInfo.zCoord - 1) * 2;  % sc_f = 2

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
    xmlwrite(fullfile(xml_folder, ['GMEMfinalResult_frame' ind '.xml']), docNode);
end