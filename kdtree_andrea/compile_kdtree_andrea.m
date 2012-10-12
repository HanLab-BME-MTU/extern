
curPath = pwd;
[kdtree_script_path, name, ext] = fileparts( mfilename('fullpath') );

cd(kdtree_script_path);

    mex kdtree_build.cpp

    mex kdtree_k_nearest_neighbors.cpp
    mex kdtree_nearest_neighbor.cpp
    
    mex kdtree_range_query.cpp    
    mex kdtree_ball_query.cpp
    
    mex kdtree_delete.cpp
    
cd(curPath);