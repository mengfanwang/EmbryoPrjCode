function display_track_from_root(movieInfo, root_id, ...
    refine_res, embryo_vid, graph_disp)
% debug use: display the linking condition of a track
%
if nargin == 4
    graph_disp = true;
end

cur_track = root_id;
pointer_loc = 1;
while pointer_loc <= length(cur_track)
    cur_parent = cur_track(pointer_loc);
    pointer_loc = pointer_loc + 1;
    cur_track = cat(1, cur_track, movieInfo.kids{cur_parent});
end
cur_track = unique(cur_track);
if length(cur_track) == 1
    fprintf('track length is 1 with node %d\n', cur_track);
    return;
end
s = [];
t = [];
w = [];
if isfield(movieInfo, 'nei') && isfield(movieInfo, 'CDist')
    for i=1:length(cur_track)
        kids = movieInfo.kids{cur_track(i)};
        s = [s, cur_track(i) + kids'*0];
        t = [t, kids'];
        nei = movieInfo.nei{cur_track(i)};
        CDist = movieInfo.CDist{cur_track(i)};
        %CDist = movieInfo.Cij{cur_track(i)};
        for j=1:length(kids)
            if isempty(CDist(nei==kids(j)))
                keyboard;
            end
            w = [w, CDist(nei==kids(j))];
        end
    end
else
    for i=1:length(cur_track)
        kids = movieInfo.kids{cur_track(i)};
        s = [s, cur_track(i) + kids'*0];
        t = [t, kids'];
    end
end
if nargin > 2 && ~isempty(refine_res)
    cnt = 0;
    uni_n = unique([s, t]);
    st_n = max(uni_n) + 1;
    for i=1:length(uni_n)
        vidx = movieInfo.voxIdx{uni_n(i)};
        frs = movieInfo.frames(uni_n(i));
        ids = refine_res{frs}(vidx);
        ids = unique(ids(ids>0));
        if length(ids) > 1
            cnt = cnt + 1;
%             rep_ns = st_n+[1:length(ids)];
%             locs = find(s==uni_n(i));
%             for j=1:length(locs)
%                 s=[s, rep_ns];
%                 t = [t, rep_ns*0+t(locs(j))];
%             end
%             loct = find(t==uni_n(i));
%             for j=1:length(loct)
%                 t=[t, rep_ns];
%                 s = [s, rep_ns*0+s(loct(j))];
%             end
%             loc = cat(2, locs, loct);
%             s(loc) = [];
%             t(loc) = [];
%             st_n = st_n+length(ids);
        end
    end
    %disp(cnt);
end

[nd_st,~,ic] = unique([s, t]);
s1 = ic(1:length(s));
t1 = ic(length(s)+1:end);
if ~isempty(w)
    if graph_disp
        G = digraph(s1,t1,w);
        figure; plot(G, 'NodeLabel',nd_st,'EdgeLabel', round(G.Edges.Weight,1));%);
    end
    blue_label_scale = 3;
else
    if graph_disp
        G = digraph(s1,t1);
        figure; plot(G, 'NodeLabel',nd_st);
    end
    blue_label_scale = 3;
end

if nargin > 3 && ~isempty(refine_res) && ~isempty(embryo_vid)
    
    [h,w,zslice] = size(refine_res{1});
    frs = movieInfo.frames(cur_track);
    unifrs = unique(frs);
    segMaxProj = uint8(zeros(h,w,3,numel(refine_res)));
    for i=1:numel(refine_res)
        tmp_g = max(embryo_vid{i},[],3);
        segMaxProj(:,:,2,i) = tmp_g;
    end
    for i=1:length(unifrs)
        voxIdx = cat(1, movieInfo.voxIdx{cur_track(frs==unifrs(i))});
        [yy,xx,~] = ind2sub([h,w,zslice], voxIdx);
        tmp_fr = false(h,w);
        tmp_fr(sub2ind([h,w], yy, xx)) = true;
        tmp_g = max(embryo_vid{unifrs(i)},[],3);
        %segMaxProj(:,:,2,unifrs(i)) = tmp_g;
        tmp_b = uint8(zeros(h,w));
        tmp_b(tmp_fr) = tmp_g(tmp_fr)*blue_label_scale;
        segMaxProj(:,:,3,unifrs(i)) = tmp_b;
    end
    
    zzshow(segMaxProj);
end