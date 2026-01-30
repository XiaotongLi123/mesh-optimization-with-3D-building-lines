function [InFace_idx] = getKnnFace(KDmodel, neighbor, face, line, KNN, interval)
    % Initialize LSM
    % Input: KD-tree model, point-face adjacency table, line segment, distance threshold
    % Output: Indices of intersecting faces
    node = getNode(line,interval);
    node_num = size(node, 1);
    face_num = size(face, 1);
    InFace_idx_ori = zeros(1024, 1);
    InFace_count = 0;
    face_label = zeros(face_num, 1);

    for n = 1:node_num
        index = knnsearch(KDmodel, node(n,:), 'K', KNN);
        for i = 1:KNN
            for j = 1:neighbor(index(1,i), 1)
                if face_label(neighbor(index(1,i), j+1), 1) ~= 1
                    InFace_count = InFace_count + 1;
                    InFace_idx_ori(InFace_count, 1) = neighbor(index(1,i), j+1);
                    face_label(neighbor(index(1,i), j+1), 1) = 1;
                end
            end
        end
    end
    InFace_idx = InFace_idx_ori(1:InFace_count, 1);
end