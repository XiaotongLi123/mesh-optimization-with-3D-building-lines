addpath(genpath('Algorithm'));

clear; clc;
%------------Read data files and preprocess------------
dataSet = 'dublin';
meshType = 'mesh_10000';
para_1 = 30;
para_2 = 0.6;
node_interval = 0.03;
dense_interval = node_interval / 20;
sample_interval = node_interval * 5;

% [point, face] = readMeshOBJ(['data\', dataSet, '\', meshType, '.obj']);
% [line_set, ~] = readLineOBJ(['data\', dataSet, '\line.obj']);

[point, face] = readMeshOBJ('mesh_dublin.obj');
[line_set, ~] = readLineOBJ('line.obj');

line_num = size(line_set, 1) / 2; % Number of original lines
pts_num = size(point, 1);         % Number of original points
face_num = size(face, 1);         % Number of original faces

para_1_set = para_1 * ones(line_num, 1);
para_2_set = para_2 * ones(line_num, 1);

error = 0;                     % Count the number of erroneous lines
error_id = zeros(line_num, 1); % Store IDs of lines that throw errors
reconstruct_line = zeros(line_num, 1); % Store IDs of successfully reconstructed lines
constrained_line_set = [];     % Record constrained edges

tic;
line_set = lineSort(line_set);

for m = 1:line_num
%     try
        m
        %------------Construct k-d tree------------
        KDmodel = KDTreeSearcher(point); % Build KD-tree
        PFneighbor = createPFneighbor(point, face); % Build point-face adjacency table
        line = line_set(2*m-1:2*m, :);

        %------------LSM initializing------------
        knn_face_idx = getKnnFace(KDmodel, PFneighbor, face, line, 3, dense_interval);
        knn_face_idx = holeFilling(knn_face_idx, face, PFneighbor);
        knn_face_idx = regionConnect(point, face, knn_face_idx, PFneighbor);
        knn_face_num = size(knn_face_idx, 1);
        
        %------------line-to-mesh-projection------------
        [region_face_idx, near_point] = IsIntersectant(point, face, knn_face_idx, line, dense_interval);
        region_face_idx = holeFilling(region_face_idx, face, PFneighbor);
        region_face_idx = regionConnect(point, face, region_face_idx, PFneighbor);
        region_face_num = size(region_face_idx, 1);
        
        %------------LSM growing------------
        l2f_dist_base = zeros(region_face_num, 1);
        for i = 1:region_face_num
            l2f_dist_base(i, 1) = line2faceDist2(point(face(region_face_idx(i,1), :)', :), line);
        end
        max_dist = getQuantile(l2f_dist_base, para_2_set(m, 1));
        if max_dist < 1e-4
            max_dist = 1e-4;
        end
        range_face_idx = getKnnFace(KDmodel, PFneighbor, face, line, 250, sample_interval);
        range_face_num = size(range_face_idx, 1);
        l2f_dist = zeros(range_face_num, 1);    
        for i = 1:range_face_num
            l2f_dist(i, 1) = line2faceDist2(point(face(range_face_idx(i,1), :)', :), line);
        end
        normal = getFacesNormal(point, face(range_face_idx, :));

        extend_region_face_idx = regionExtend(face, range_face_idx, region_face_idx, normal, l2f_dist, PFneighbor, para_1_set(m, 1), max_dist);
        extend_region_face_idx = shapeRepair(point, extend_region_face_idx, face, PFneighbor);
        extend_region_face_idx = holeFilling(extend_region_face_idx, face, PFneighbor);
        extend_region_face_num = size(extend_region_face_idx, 1);

        %------------LSM clipping------------     
        [segment_region_face_idx, this_error] = regionSegment(point, face, extend_region_face_idx, region_face_idx, line, constrained_line_set, node_interval);
        if this_error == 1
            error = error + 1;
            error_id(error, 1) = m;
            continue;
        end
        [segment_region_face_idx, error] = holeFilling(segment_region_face_idx, face, PFneighbor);
        if haveConstrainedEdge(face, segment_region_face_idx, constrained_line_set)
            error = error + 1;
            error_id(error, 1) = m;
            continue;
        end
        segment_region_face_num = size(segment_region_face_idx, 1);

        %------------LSM retriangulation------------
        [face, point, normal_angle, constrained_line_set, ~, this_error] = reconstruct_ours(face, point, segment_region_face_idx, line, PFneighbor, constrained_line_set, node_interval); % Reconstruct to get new face table
        new_face_num = size(face, 1);

        if this_error == 1
            error = error + 1;
            error_id(error, 1) = m;
            continue;
        else
            face_num = new_face_num;
            pts_num = size(point, 1);
            reconstruct_line(m, 1) = 1;
        end
%     catch
%         error = error + 1;
%         error_id(error, 1) = m;
%         continue;
%     end
end
toc;
writeMeshOBJ(['data\', dataSet, '\', meshType, '_', num2str(para_1), '_', num2str(para_2), '.obj'], point, face);

function [result] = getQuantile(data, q)
    data_sort = sort(data, 'ascend');
    idx = ceil(size(data, 1) * q);
    result = data_sort(idx, 1);
end