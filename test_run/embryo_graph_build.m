file_path = '/home/congchao/Desktop/mats';

detections = cumulateXMLs(file_path);

[dat_in, excess_node, c_en, c_ex] = build_embryo_graph(detections);
save('embryo_10tb_heuristic_graph.mat', 'dat_in', 'excess_node','c_en', 'c_ex', '-v7.3');






