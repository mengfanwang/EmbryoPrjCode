function movieInfo = removeCellRelation(movieInfo, root_cell_id, nei_cell_ids)
% remove the neighborhood relationship among root_cell_id and other cells
% included in nei_cell_ids

f0 = movieInfo.frames(root_cell_id);
for i=1:length(nei_cell_ids)
    if nei_cell_ids(i) > 0 && nei_cell_ids(i) <= numel(movieInfo.frames)
        f1 = movieInfo.frames(nei_cell_ids(i));
        if f1>f0
            od = find(movieInfo.nei{root_cell_id}==nei_cell_ids(i));
            if ~isempty(od)
                if length(movieInfo.nei{root_cell_id}) == 1
                    movieInfo.nei{root_cell_id} = [];
                    movieInfo.CDist{root_cell_id} = [];
                    movieInfo.CDist_i2j{root_cell_id} = [];
                    movieInfo.Cij{root_cell_id} = [];
                    movieInfo.ovSize{root_cell_id} = [];
                else
                    movieInfo.nei{root_cell_id}(od) = [];
                    movieInfo.CDist{root_cell_id}(od) = [];
                    movieInfo.CDist_i2j{root_cell_id}(od,:) = [];
                    movieInfo.Cij{root_cell_id}(od) = [];
                    movieInfo.ovSize{root_cell_id}(od) = [];
                end
                
                od = find(movieInfo.preNei{nei_cell_ids(i)}==root_cell_id);
                if length(movieInfo.preNei{nei_cell_ids(i)}) == 1
                    movieInfo.preNei{nei_cell_ids(i)} = [];
                    movieInfo.Cji{nei_cell_ids(i)} = [];
                    movieInfo.CDist_j2i{nei_cell_ids(i)} = [];
                    movieInfo.preOvSize{nei_cell_ids(i)} = [];
                else
                    movieInfo.preNei{nei_cell_ids(i)}(od) = [];
                    movieInfo.Cji{nei_cell_ids(i)}(od) = [];
                    movieInfo.CDist_j2i{nei_cell_ids(i)}(od,:) = [];
                    movieInfo.preOvSize{nei_cell_ids(i)}(od) = [];
                end
            end
        elseif f1<f0
            od = find(movieInfo.nei{nei_cell_ids(i)}==root_cell_id);
            if ~isempty(od)
                if length(movieInfo.nei{nei_cell_ids(i)}) == 1
                    movieInfo.nei{nei_cell_ids(i)} = [];
                    movieInfo.CDist{nei_cell_ids(i)} = [];
                    movieInfo.CDist_i2j{nei_cell_ids(i)} = [];
                    movieInfo.Cij{nei_cell_ids(i)} = [];
                    movieInfo.ovSize{nei_cell_ids(i)} = [];
                else
                    movieInfo.nei{nei_cell_ids(i)}(od) = [];
                    movieInfo.CDist{nei_cell_ids(i)}(od) = [];
                    movieInfo.CDist_i2j{nei_cell_ids(i)}(od,:) = [];
                    movieInfo.Cij{nei_cell_ids(i)}(od) = [];
                    movieInfo.ovSize{nei_cell_ids(i)}(od) = [];
                end
                od = find(movieInfo.preNei{root_cell_id}==nei_cell_ids(i));
                if length(movieInfo.preNei{root_cell_id}) == 1
                    movieInfo.preNei{root_cell_id} = [];
                    movieInfo.Cji{root_cell_id} = [];
                    movieInfo.CDist_j2i{root_cell_id} = [];
                    movieInfo.preOvSize{root_cell_id} = [];
                else
                    movieInfo.preNei{root_cell_id}(od) = [];
                    movieInfo.Cji{root_cell_id}(od) = [];
                    movieInfo.CDist_j2i{root_cell_id}(od,:) = [];
                    movieInfo.preOvSize{root_cell_id}(od) = [];
                end
            end
        end
    end
end

end