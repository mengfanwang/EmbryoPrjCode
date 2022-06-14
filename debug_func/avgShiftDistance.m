function [shifts, avgShifts] = avgShiftDistance(movieInfo)
shifts = cell(max(movieInfo.frames), 1);
%avgShifts = zeros(max(movieInfo.frames), 1);
for i=1:numel(movieInfo.tracks)
    for j=2:length(movieInfo.tracks{i})
        n1 = movieInfo.tracks{i}(j-1);
        n2 = movieInfo.tracks{i}(j);
        f1 = movieInfo.frames(n1);
        f2 = movieInfo.frames(n2);
        if f2-f1==1
            od = movieInfo.nei{n1}==n2;
            shifts{f2} = cat(1, shifts{f2}, movieInfo.CDist{n1}(od));
        end
    end
end

avgShifts = cellfun(@mean, shifts);
end
