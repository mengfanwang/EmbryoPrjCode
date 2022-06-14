Main function: track_refine_YinanData.m

Part 1: Initialization
	Parameter:	embryo_vid_org: data
				org_eigMaps: {1} 2d curvature {2} 3d curvature
				sc_f: scale factor of downsampling

	Function 1.1 graphPara_cell: initialize parameters for tracking
		g:
		.particleNum ----------	number of cells
		.use_translation_as_drift whether use translation as drift
		.translation ----------	translation drift
		.maxIter --------------	maximum iteration (of what?)
		.realEnter ------------ chi2inv(1-0.01/g.particleNum, 1)/2. Why?
		.c_en/c_ex ------------ enter and exit cost
		.observationCost ------ observation cost
		.k -------------------- maximum number of jump allowed


Part 2: Tracking
	|---2.1	mcfTracking_cell
		|---2.2 InitialDriftEstimate
		|---2.3 tree2tracks_cell
			|---2.4 detection2tracks: not used？
		|---2.5 transitCostInitial_cell
			|---2.6 ovDistanceMap

	movieInfo:
		.x/y/zCoord ------- xyz locations
		.n_perframe -------	number of ids per frame
		.label ------------	id list of all frames
		.frames ----------- t locations
		.vox --------------	voxel location of all detections
		.voxIdx ----------- voxel index of all detections
		.tracks -----------
		.particle2track ---
		.orgCoord --------- original xyz coordinate
		.Ci --------------- observation cost

	Function 2.1 mcfTracking_cell: tracking based on min-cost flow/circulation
	Input:	det_maps ----------	segmentation id maps
			embryo_vid -------- data
			threshold_res -----	threshold of segmentation
			varMaps/eigMaps/g/q
	Output: movieInfo ---------
			movieInfoAll ------
			refine_res --------
			refine_resAll -----
			threshold_res -----
			threshold_resAll --
	Note: 1.the size of _All is maxIter*3+2. Why?


	Function 2.2 InitialDriftEstimate: estimate the drift by aligning the data with translation
	Input: 	time_ordered_data - data
	Output:	drift

	Function 2.3 tree2tracks_cell: build initial tracking info 
	Input:	det_maps ----------	segmentation id maps
			detections -------- default as false
			q
	Output:	detections_YXZ ---- (yxzt) center location of detections
			voxIdxList -------- voxel index of all detections
			movieInfo:
			.x/y/zCoord ------- xyz locations
			.n_perframe -------	number of ids per frame
			.label ------------	id list of all frames
			.frames ----------- t locations
			.vox --------------	voxel location of all detections
			.voxIdx ----------- voxel index of all detections
			.tracks ----------- not used?
			.particle2track ---	not used?
	Note: 1.vox will be influenced by drifting but voxIdx will not. This is to make the location identificatoin easier

	Function 2.5 transitCostInitial_cell: initial variance and mean of distances between observed object locations and estimated locations from data
	Input: 	det_maps voxIdxList movieInfo g
			dets_YXZ ----------	(yxzt) center location of detections
	Output: movieInfo

















------------------------------ Element functions --------------------------------------------
data_scaling: (crop and) scale data
Input:	sc_f ------------------	scale factor
		st_loc ----------------	start location
		sz_crop ---------------	size of crop data
		org_refine_res -------- segment result
		embryo_vid_org -------- data
		gt_mat_org ------------ ?? (may not be used)
		org_threshold_res ----- segment threshold
		org_varMap ------------ variance map
		org_eigMaps ----------- curvature map
Output	refine_res/embryo_vid/gt_mat/threshold_res/varMap/eigMaps