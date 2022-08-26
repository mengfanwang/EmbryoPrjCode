Main function: track_refine_YinanData.m

Part 1: Initialization
	Parameter:	embryo_vid_org: data
				org_eigMaps: {1} 2d curvature {2} 3d curvature
				sc_f: scale factor of downsampling

	Function 1.1 graphPara_cell: initialize parameters for tracking
		g:
		.particleNum ----------	total number of cells
		.use_translation_as_drift whether use translation as drift (Function 2.2)
		.translation ----------	translation drift
		.maxIter --------------	maximum iteration (of what?)
		.realEnter ------------ chi2inv(1-0.01/g.particleNum, 1)/2. Why?
		.c_en/c_ex ------------ enter and exit cost (= realEnter)
		.observationCost ------ observation cost
		.k -------------------- maximum number of frames of jump allowed
		.cycle_track ---------- = true using mcc = false using mcf
		.n_nodes -------------- nodes of graph = 2*particleNum + 1/2
		.excess_node ---------- [1 g.nodes] source and sink for mcf graph
		.trackLength4var ------ tracks with smaller length will not be used for estimation
		.jumpCost ------------- moveInfo.jumpRatio
		.applyDrift2allCoordinate =True update coordinate location = False update cost in 2.19
		.useOldArcProps ------- use old movieInfo.arc_avg_mid_std in 2.16

Part 2: Tracking
	|---2.1	mcfTracking_cell
		|---2.2 InitialDriftEstimate
		|---2.3 tree2tracks_cell
			|---2.4 detection2tracks
		|---2.5 transitCostInitial_cell
			|---2.6 ovDistanceMap
			    |---2.7 verifyInData
			    |---2.8 ovDistanceRegion
		    |---2.9 fitTruncGamma
		    |---2.10 overlap2cost
	    |---2.11 trackGraphBuilder_cell
	    |---2.12 mccTracker
	    	|---2.13 mcc4mot
	    	|---2.14 mcc4mot_cell
	    	|---2.15 upt_cost_with_Dist
	    	|---2.16 stable_arc_cost_extract
	    		|---2.17 arc_cost_in_track 
	    	|---2.xx relinkJumpOveredCell
	    |---2.xx transitCostUpt_cell
	    	|---2.xx driftFromTracks
	    	|---2.xx getStdFromTracks_cell   -->2.9




	movieInfo:
		.x/y/zCoord ------- xyz locations
		.n_perframe -------	number of ids per frame
		.label ------------	id list of all frames
		.frames ----------- t locations
		.vox --------------	voxel location of all detections
		.voxIdx ----------- voxel index of all detections
		.orgCoord --------- original xyz coordinate
		.Ci --------------- Initial: observation cost

		.ovGamma ---------- parameters of estimated Gamma distribution
		.CDist ------------ max overlap distance
		.CDist_i2j -------- both two overlap distances
		.CDist_j2i -------- backward distance (the latter of i2j)
		.Cij -------------- forward  cost (distance -> z-score)
		.Cji -------------- backward cost
     	.nei -------------- forward  neighbor indexes
     	.preNei ----------- backward neighbor indexes
     	.ovSize ----------- forward  overlap size
     	.preOvSize -------- backward overlap size 	

		.parents ---------- parents of each detection
		.kids ------------- kids of each children
		.track_bf_merge --- trajectories before merging tracks




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


	Function 2.2 InitialDriftEstimate: estimate the drift by imregtform
	Input: 	data
	Output:	drift

	Function 2.3 tree2tracks_cell: build initial tracking info 
	Input:	det_maps ----------	segmentation id maps
			detections -------- 
			q
	Output:	detections_YXZ ---- (yxzt n*4) center location of detections 
			voxIdxList -------- voxel index of all detections
			movieInfo:
			.x/y/zCoord ------- xyz locations
			.n_perframe -------	number of ids per frame
			.label ------------	id list of all frames
			.frames ----------- t locations
			.vox --------------	voxel location of all detections
			.voxIdx ----------- voxel index of all detections
	Note: 1.vox will be influenced by drifting but voxIdx will not. This is to make the location identificatoin easier (forget the reason)

	Function 2.5 transitCostInitial_cell: initial variance and mean of distances between observed object locations and estimated locations from data
	Input: 	det_maps voxIdxList movieInfo g
			dets_YXZ ----------	(yxzt n*4) center location of detections 
	Output: movieInfo
			.ovGamma ---------- parameters of estimated Gamma distribution
			.CDist ------------ max overlap distance
			.CDist_i2j -------- both two overlap distances
			.CDist_j2i -------- backward distance (the latter of i2j)
			.Cij -------------- forward  cost (distance -> z-score)
			.Cji -------------- backward cost
	     	.nei -------------- forward  neighbor indexes
	     	.preNei ----------- backward neighbor indexes
	     	.ovSize ----------- forward  overlap size
	     	.preOvSize -------- backward overlap size 

	Function 2.6 ovDistanceMap: calculate the distance between regions across frames based on overlapping ratio
	Input:  curMap ------------ id map of current frame
			nextMap ----------- id map of future frame 
			curVoxIdxCells ---- voxel indexes of current frame
			nextVoxIdxCells --- voxel indexes of future frame
			drift ------------- translation between two frames
	Output: neighbors --------- cells overlapped in the future frame
			distances --------- max average distance between two regions
			distances_dirwise - average min distance between two regions
	Note: "future" rather than "next" because of jump

	Function 2.7 verifyInData: remove locations out of the boundary.
	Input:  yxz --------------- locations
	        sz_yxz ------------ boundary
    Output: out --------------- locations in the boundary 

    Function 2.8 ovDistanceRegion: cal the overlapping distance between two regions
    Logic:  if ovflag
    		    distances = -log IOU
    		else
                average min distance discussed in Dissertation P75.
    Input:  curRegVox --------- current region voxel locations
    		nextRegVox -------- future region voxel locations
    		frame_shift ------- shift
    		ovFlag ------------ whether using overlapping ratio as distance. default = false
	Output: maxDistance ------- the larger distance 
			minDistance ------- the smaller distance
			distances --------- average min distance between two regions
			re_ratio ---------- (boxSize - future region) / boxSzie

	Function 2.9 fitTruncGamma: fit Truncated Gamma distribution
	Input:  data -------------- the nearest distance
	Output: phat -------------- Gamma parameter
	Note: In 2.8, distance of two small regions is 100. It may influence the fitting.

	Function 2.10 overlap2cost: change the overlap distance to z-score (cost)
    Input:  ov_score ---------- input score (distance)
            phat -------------- Gamma parameter
            jump_punish ------- jump ratio. jumps are punished.
    Output: cost -------------- converted z-score
    		oc ---------------- sqrt(cost)
	Note: probability * jump is weird. How about (p)^(1/punish) ?

    Function 2.11 trackGraphBuilder_cell: build tracking graph
    Input:  movieInfo g
    Output: g
			orgG -------------- work on mcf graph. digraph of data_in
            dat_in ------------ [s t w] for mcf and {detection, transition} for mcc

    Function 2.12 mccTracker: get trajectories by solving min-cost circulation
    Input:  dat_in ------------ {detection, transition, (append)}_arcs
            movieInfo g 
    Output: movieInfo
			.parents ---------- parents of each detection
			.kids ------------- kids of each children
			.track_bf_merge --- trajectories before merging tracks
			.tracks ----------- trajectories
			.pathCost --------- nan
			.particle2track --- [track it belongs to, position in the track, track_bf_merge it belongs to]
			.jumpRatio -------- ratio of different distances of jumps

   	Function 2.13 mcc4mot: solving original min-cost circulation
   	Input:  detection_arcs,transition_arcs
   	Output: trajectories, costs

   	Function 2.14 mcc4mot_cell: solving mcc with appended arcs(split/merge)
   	Input:  detection_arcs,transition_arcs, append_arcs
   	Output: trajectories, costs
   			track_bf_merge ---- trajectories before merging tracks
   			parents ----------- parent nodes of the current nodes
   			kids -------------- kid nodes of the current nodes
	Note: need double-check (line 158-181)

	Function 2.15 upt_cost_with_Dist: update cost if something changed
	Input:  movieInfo.nei/ovGamma/CDist/jumpRatio
	Output: movieInfo.Cij/Cji

	Function 2.16 stable_arc_cost_extract:
	Input:  movieInfo, g
	Output: movieInfo
	        .arc_avg_mid_std --	[mean, median, std, location] of a track


    Function 2.17 arc_cost_in_track: extract the costs from a track
    Input:  movieInfo, track_id
    		split_flag --------	not implemented and used
    		max_jump ----------	allowed maximum jumps


	Function 2.18 transitCostUpt_cell: update transition cost using existing trajectories
	Input:  movieInfo, g
	Output: movieInfo

	Function 2.19 driftFromTracks: use common motion of adjacent frames to estimate drift
	Input:  movieInfo, g
	Output: dfEst -------------	estimated drift
			movieInfo.drift --- dfEst

	Function 2.20 getStdFromTracks_cell: update the gamma distribution parameters
	Input: 	movieInfo, g
	Output: phatDist ---------- gamma parameters
			










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