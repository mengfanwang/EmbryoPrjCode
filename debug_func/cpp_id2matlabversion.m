function res = cpp_id2matlabversion(movieInfo, idx)

res = idx + 1;
cumulate = 2*movieInfo.n_perframe;
for i = 2:length(cumulate)
    cumulate(i) = cumulate(i) + cumulate(i-1);
end
for i = 1:length(res)
    f = find(cumulate>res(i), 1) - 1;
    if f > 1
        res(i) = res(i) - cumulate(f) / 2;
    end
end





    