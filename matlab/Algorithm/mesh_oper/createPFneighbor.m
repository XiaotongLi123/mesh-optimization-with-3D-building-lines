function [PFneighbor] = createPFneighbor(point, face)
    % Construct point-face adjacency table
    % Input: 3D points (either point IDs or coordinates, mainly used to get the number of points; C++ version may consider overloading), faces
    % Output: Face adjacency list for each point

    pts_num = size(point, 1);
    face_num = size(face, 1);
    
    % Initialize adjacency table
    % Each row: first element stores the number of adjacent faces (max 19)
    % Remaining elements store the indices of adjacent faces, padded with 0 if needed
    PFneighbor = zeros(pts_num, 20); 
    
    for i = 1:face_num
        for j = 1:3
            point_face_neigh = PFneighbor(face(i,j), 1); % Current number of adjacent faces for this point
            if point_face_neigh >= 19
                continue; 
            end
            PFneighbor(face(i,j), point_face_neigh + 2) = i; % Add adjacency relation at the corresponding position
            PFneighbor(face(i,j), 1) = PFneighbor(face(i,j), 1) + 1; % Update adjacency count
        end
    end
end