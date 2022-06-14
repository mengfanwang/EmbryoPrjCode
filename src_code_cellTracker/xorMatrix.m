function [out_m1, out_m2] = xorMatrix(in_m1, in_m2)
% remove the rows in two matrics with overlapped element

vec_m1 = in_m1(:);
vec_m2 = in_m2(:);
if isempty(intersect(vec_m1, vec_m2))
    out_m1=in_m1;
    out_m2=in_m2;
    return;
end
flag = true(size(in_m1,1), 1);
for i=1:size(in_m1,1)
    if ~isempty(intersect(in_m1(i,:),vec_m2))
        flag(i) = false;
    end
end
out_m1=in_m1(flag,:);

flag = true(size(in_m2,1), 1);
for i=1:size(in_m2,1)
    if ~isempty(intersect(in_m2(i,:),vec_m1))
        flag(i) = false;
    end
end
out_m2=in_m2(flag,:);

end