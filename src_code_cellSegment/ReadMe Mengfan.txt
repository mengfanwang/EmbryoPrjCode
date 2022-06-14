Main function: cellSegment4YinanData.m

Part 1: Initialization

Part 2: Use SynQuant to get a preliminary result
	Paramter: sigma minIntensity

	Function 2.1 Synquant4Embryo_Paramater: get detection result by SynQuant
	Input:  sm_im ---------- smoothed data
			q -------------- q.minIntensity
	Output: zMap ----------- z-score map of detections
			synId ---------- id map of detections
			fMap ----------- region > minIntensity
	Note: synId(zMap <=1 | ~fMap) = 0 

	Output: z_mat ---------- cell of zMap
			id_mat --------- cell of synId
			fMaps ---------- cell of fMap

Part 3: Calculate 2D&3D principle curvature
	Parameter: sigma (different for 2D&3D)

	Function 3.1 principalCv2d(3d): calculate principal curvature
	Input:  vid ------------ data
			synId ---------- id map of detections
			sigma ---------- smoothness (length = 1 for 2d & 3 for 3d)
	Output: eig_all -------- maximum eigenvalue of Hessian matrix after smooth
			overlay_cl ----- visualization of eigenvalues larger than valid_eig over on synId

	Output: eig_res_2d ----- cell of 2d eigenvalue map
			eig_res_3d ----- cell of 3d eigenvalue map	
			eig_overlay ---- cell of 3d overlay_cl


Part 4: Calculate the noise variance map
	Parameter: scale_term

	Function 4.1 calVarianceStablizationBY: estimate noise variance
	Input:  imregion1 -------- data scaled to around [0,255]
		    varRatio --------- the percentage of selecting sigma2
		    gap -------------- gap between target pixel and surrounded pixels
    Output:
    	    varMap ----------- pixel-wise noise variance map
    	    sigma2 ----------- variance at varRatio
    	    xx --------------- variance curve	

    Output: varMap ----------- 3*2 cell:
							  Before stabilization     After stabilization
							         varMap                  varMap 
							         sigma2                  sigma2
							         xx                      xx
    Note:
    1. Data is scaled to around [0,255] by the scale term 
    2. The estimation after stabilization is quite inaccurate on registered data which don't follow Possion-Gaussian distribution. Use ConvexVST as a substitution.


Part 5: Refine detection results combining curvature and noise variance information
	|---5.1 regionWiseAnalysis4d
		|---5.2 initial_q
		|---5.3 initial_Orst
		|---5.4 region_sanity_check
			|--- 5.5 rearrange_id
		|---5.6 refineOneRegion_with_seed
		 	|---5.7 get_local_area
		 	  	|---5.8 cropNeedRegions
		 	  	   	|--- crop3D
		 	  	   	|--- 5.9 pickLargestReg
		 	  	   	|--- 5.10 coordinate_transfer
		 	|---5.11 fgDetectSynQuant
		 	  	|---5.12 ordstat4fg
		 	  		|---5.13 regionSig
		 	  			|---5.14 findValidNeiReg
		 	  			|---5.15 ordStats
		 	  				|---5.16 ordStatApproxKsecWith0s_mat
 	  				|---5.17 refine_with_seedRegion
		 	  	|--- scale_image
		 	  	|---5.18 splitFGintoCells
		 	  		|---5.19 regionGrow
		 	  			|---5.20 graphCut_negHandle_mat
		 	  				|---neighbours_mat
		 	  		|---5.21 extraVoxReassign
		 	  		|---5.22 seed_map_gen
		 	|---5.23 segmentCurrentRegion
		 	  	|---5.24 gapTest3dV2
		 	  	   	|---5.19 regionGrow
		 	  	   	|---5.25 edgeTest3dV2
		 	  	   		|---5.26 gapRefine
		 	  	   		|---5.27 neighbor_pixels
		 	  	   		|---ordStatApproxKsec
		 	  	   		|---ordStatApproxKsecWith0s
		 	|---5.28 region_refineV2
		 		|---5.19 regionGrow
		 		|---5.21 extraVoxReassign
		 		|---5.29 shrinkTest3d
		 	|---5.30 fgDetectSynQuant_thresGiven

	
	Function 5.1 regionWiseAnalysis4d:
	Input:  idMap ----------- id map of detections
		    eigMap ---------- {2d eigenvalue map, 2d eigenvalue map}
		    vid ------------- data scaled to around [0,255]
		    varMap ---------- variance map
		    test_ids --------
		    tif_id ---------- for write image
    Output: newIdMap --------
    		thresholdMap ----

    Function 5.2 initial_q: initialize parameters for region growing
    	q:
    	.minSeedSize -------- objects smaller than it are removed in region_sanity_check
    	.shift -------------- shift size for cropping local area
    	.fgBoundaryHandle --- method to process fg/boundary touch 'leaveAloneFirst' 'repeat' 'compete'
    	.cost_design -------- two elements for designing the edge cost of region growing graph
    						  1.method: = 1 use the average of two voxels 2/(p1+p2)^exponent 
    						  			= 2 use sqrt root of two voxels 1/sqrt(p1*p2)^exponent
    						  2.exponent						
	  	.neiMap ------------- connectivity for segmentCurrentRegion
	  	.growConnectInRefine  connectivity for gapTest3dV2


	Function 5.3 initial_Orst: initialize parameters for fg detection and gap significance test
		Orst:
		.NoStbVarMap -------- variance map before stabilization
		.NoStbVar ----------- estimated variance before stabilization
		.stbVarMap ---------- variance map after stabilization
		.imProcMethod ------- 'stb' or 'noStb': use stabilization data or original data
		.corrTerm ----------- = max(stb)/255, stabilized data scale term
		.fgTestWay ---------- foreground hypothesis testing method (order-statistics or others)
		.gapTestWay --------- gap testing method (local ost or others)
		.refine_cc ---------- use the region containing seed only or not 
		.fgVarSamplewiseEst - use pixel-wise variance map 
		.curStbVar ---------- average variance of current region

	Function 5.4 region_sanity_check: remove objects smaller than minSz
	Input:  regLabel -------- id map of detections
			minSz ----------- removing threshold
	Output: regLabel -------- id map after removing and sorting
			redundant_flag -- = True if removed

	Function 5.5 rearrange_id

	Function 5.6 refineOneRegion_with_seed:
	Input:  seed_id --------- id of current region
			yxz ------------- index list of current region
			vid ------------- data
			vid_stb --------- data after stabilization
			idMap ----------- id map of detections
			eig2d/3d -------- 2d/3d eigenvalue maps
			Orst/q ---------- parameters
	Output: newLabel -------- refined label map
			comMaps --------- output of 5.8 cropNeedRegions. information about the local region
			fgReDo ---------- =True means the foreground is too small and redo the pipeline

	Function 5.7 get_local_area: 
	Input: the same as 5.6 (i is seed_id)
	Output: comMaps --------- output of cropNeedRegions
 	        OrSt:
 	        .stbVarCropMap -- cropped .stbVarMap
 	        .NoStbVarCropMap  cropped .NoStbVarMap

    Function 5.8 cropNeedRegions: crop a region defined by yxz and shift
    Parameter: sm_term = [1 1 1] for vid_sm
    Input:  vid vid_stb idMap eig2d/3d yxz: the same as 5.6 refineOneRegion_with_seed
    		reg_id ---------- seed_id
    		improc ---------- processing method ('stb' or 'noStb')
    		shift ----------- shift size for cropping the region
    Output: comMaps
    		.vidComp -------- cropped vid
    		.linerInd ------- index of the cropped region in the image 
    		.yxz_edge_Flag -- = True if touching boundary of the image
    		.vidStbComp ----- cropped vid_stb
    		.idComp --------- cropped idMap
    		.vid_sm --------- smoothed (stabilized) data depend by improc
    		.regComp -------- the largest component of idComp (seed region)
    		.fmapCompInit --- candidate foreground map
    		.fmapComp ------- foreground map after post-processing
    		.eigPosMap ------ (binary) positive 3d eigenvalue map
    		.score3dMap ----- positive 3d eigenvalues
    		.eig2dPosMap ---- (binary) positive 2d eigenvalue map
    		.score2dMap ----- positive 2d eigenvalues
    		.newIdComp ------ cropped id map updated in post-processing

	Function 5.9 pickLargestReg: pick up the largest connected component
	Input:  bw -------------- binary image
			connect --------- connectivity
			valid_idx ------- (option) only these indexes are considered
	Output: out_bw ---------- binary image containing the largest component
			rm_flag --------- = True if removing anything
			idx ------------- index of the component

	Function 5.10 coordinate_transfer: change the in_idx from orgSz to cropSz
	Input: 	in_idx, orgSz, cropStart, cropSz
	Output: out_idx

	Function 5.11 fdDetectSynquant: detect the foreground within a cropped region, and growing to the appropriate boundary by graph cut
	Input: 	comMapsIn ------- output of 5.8 cropNeedRegions
			OrSt/q ---------- parameters
	OutPut:	comMaps
			.pickedThreshold  threshold to get foreground
			.score2d/3dmap -- scale so that fg (from ordstat4fg) part is [1e-3,1]
			.fmapCompInitBndIdx xy-axes boundary of fmapCompInit
			.repeatedSeed --- the seed region use for 'repeat' option, = fg
			.p_thres -------- significance threshold for gap test

	Function 5.12 ordstat4fg: find the most significant region in the given region
	Parameter: suplb: 30/255 of intensity; infub: 5/255 of intensity
	Input:	vidComp --------- data
			vid_sm ---------- smoothed data
			seedRegion ------ seed region
			other_id_map ---- id map but not seed region
			fmapComp -------- candidate foreground map
			OrSt ------------ parameters
			minSize --------- minimal size						
	Output: fg -------------- foreground
	        threshold ------- threshold to get the foreground
    Note:
   	1. other id map is dilated and removed (Line 54). Is it possible to remove the result? It may also contain some part of foregrounds.

   	Function 5.13 regionSig: select fg/bg and then computer order-statistics score
   	Parameter:  bnd_nei_radius = 3: radius of background
   	Input:  regMap ---------- foreground above threshold
   			vidMap ---------- data
   			fmap ------------ candidate foreground map
   			validNeiMap ----- candidate background
   			OrSt ------------ parameters
	Output: zscore/sigma 
			z_debias -------- = 0 if use ost
			sz_nei ---------- size of background
	Note: compare regMap and validNeiMap, both are in fmap


	Function 5.14 findValidNeiReg: find the neighboring area for testing significance of foreground 
	Input: 	regMap ----------- foreground above threshold
			validNeiMap ------ candidate background
			gapSz ------------ the num of cycled neighbors that not used for neighboring detection
			maxRadius -------- the maximum radius background can have
			radiusBasedFlag -- true means the neighboring area has the same radius as foreground region
			se_max_radius ---- strel with the maximum radius
	Output: nei_loc ---------- candidate background (neighborhood)
	Note:
	1. radius of background = min(sqrt(num)/pi, maxRadius). When radiusBasedFlag = True, nei_loc = imdilate(fgMap, se) & validNeiMap; Otherwise, dilate the background until the size is larger than foreground.


	Function 5.15 ordStats: computer the order-statistics score with different ways
	Input: 	fg --------------- foreground
			fg_neighbors ----- background
			nanVec ----------- other parts
			Orst ------------- parameters
	Output: zscore/sigma
			z_debias --------- = 0 if use ost
	Note: varMap is not used when using ost


	Function 5.16 ordStatApproxKsecWith0s_mat: input contain three parts, fg/bg and others. Assume the three parts form a distribution, but consider the difference of fg/bg only


	Function 5.17 refine_with_seedRegion: refine the foreground (1) keep the largest part (2) erode from z-axis (3) remove objects small than minSize
	Input:	fg seedRegion minSize
	Output:	fg

	Function 5.18 splitFGintoCells: give fg from order statistics, we test if the fg contains several seeds and if so, split fg based on the seeds
	Input:	fg ---------------- foreground given by ost
			fgboundary --------	boundary of fmap
			seedRegionIn ------ seedRegion
			comMaps 
			append_id --------- id of other regions don't change
			q
	Output:	comMaps
			.newIdComp -------- new seed map from seed_map_gen

	Function 5.19 regionGrow: grow regions to their boundary based on graph-cut
	Input:	newLabel ----------	label of regions	
			scoreComp ---------	score map, principal curvature of gradient
			fMap -------------- foreground can grow
			connect -----------	connectivity
			cost_design ------- from q
			bg2sink -----------	= True means connect background to sink
	Output: newLabel ----------	new label after growing
			n ----------------- max label

	Function 5.20 graphCut_negHandle_mat: build a graph for region growing
	Input: 	vid ---------------	score map, principal curvature of gradient
			fMap --------------	foreground can grow
			s/TMap ------------	voxels belong to s/t
			connect ----------- connectivity (6/26)
			cost_design -------	from q			
			bg2sink -----------	= True means connect background to sink
	Output: dat_in ------------	[out in cost] edge list
			src/sink_node ----- node id of source/sink

	Function 5.21 extraVoxReassign: reassign voxels not classified when growing
	Input:	newLabel ---------- labels of regions
			fg ---------------- foreground
	Output: newLabel
	Note: reassign on 2d firstly and then 3d

	Function 5.22 seed_map_gen: generate the seed maps from existing foreground based on curvature map
	Input:	fg ---------------- foreground
			gap3d/2d ---------- 3d/2d curvature map
			min_seed_sz ------- min seed size	
	Output: seedsMap ---------- seed map
	Note: remove 3d positive curvatures firstly. If all removed, just keep because the region is very small or not cell. Otherwise, refine the result based on 2d curvature map.
	
	Function 5.23 segmentCurrentRegion: test the gaps in the region to see if it can be segmented
	Input:	comMaps reg_id q OrSt
	Output:	newLabel edgeLabel  output of 5.25 edgeTest3dV2
			testFlag ----------	flag of running gap test (5.24) or not

	Function 5.24 gapTest3dV2: test if the gap is significant enough in a 3D manner
	Input:	L ----------------- label of seed map
			comMaps reg_id q OrSt
	Output:	newLabel regTest   	output of 5.25 edgeTest3dV2

	Function 5.25 edgeTest3dV2: test if the edge is indeed weaker than surrounding area or simply because of noise
	Input:	vid --------------- original 3d image data
			label ------------- current segmentation results
			fgMap ------------- current foreground map, which contains all boundary area whose principle
								curvature is also large
			connect ----------- the way to determine neighbors
			p_thres ----------- the p-value threshold to define significant edge
			OrSt
	OUTPUT: newLabel ---------- new segmentation which merges some segments by removing all fake edges 			regComp -----------	binary data labeling gap

	Function 5.26 gapRefine: there two triangles in 5.25 edgeTest3dV2 because of the method to detect gap. This step is to remove the triangles by detecting 4 furthest points in the cell&gap overlap.
	Input: 	slabel ------------ label map
			idpair ------------	the pair of id with gap
			edge_locs ---------	the locations of pixels in gap
	Output: edge_locs

	Function 5.28 region_refineV2: refine the region label with shrink and graph-cut
	Input:	inLabel ----------- label map
			comMaps q
			bg2sinkLink ------- link bg to sink in gap detection. default false
	Output:	newLabel ---------- label map

	Function 5.29 shrinkTest3d: if its real cell, not matter how we shrink, it should still be a unified
	Input: 	inLabel ------------ label map
			shrinkScale -------- the scale to shrink the data
			comMaps q
	Output:	outLabel ----------- label 
			reg_split ---------- = True if the region is split after shrink

	Function 5.30 fgDetectSynQuant_thresGiven: almost same as 5.11 but using given threshold rather than testing


------------------------------  Mengfan's code   --------------------------------------------
foregroundRegistration: register foreground, label map and z_score map with given transformation matrix


------------------------------ Element functions --------------------------------------------
crop3D: crop a region containing objects from 3D image
Input:  vidIn ------------------ the image to crop from
        yxz -------------------- index or subscript of objects 
        shift ------------------ shift size for cropping the region
Output: vidOut ----------------- the cropped region
		vIdx ------------------- index of the cropped region in the image 
		yxz_edge_Flag ---------- = True if touching boundary of the image
		loc_org_xyz ------------ = [xmin, ymin, zmin], the lower bound of the cropped region

scale_image(im, vlow, vhigh, ilow, ihigh): scale the input image IM to the range [vlow, vhigh] from the range [ilow, ihigh]
imo = (im-ilow)/(ihigh-ilow) * (vhigh-vlow) + vlow

ind2sub/sub2ind`_direct: the same but supporting 3D matrix directly

neighbours_mat: find the neighbor of given indexes
Input: 	ind ---------------	given indexes
		YourMatrix --------	given matrices
		connect ----------- connectivity
Output: nei_mat ----------- the neighbor indexes
Note: nei_mat size is size(ind)*connect, neighbors out of the boundary is NaN.

Function 5.5 rearrange_id: reorder id map to keep continuity
Input: 	idMap
Output: idMap
		cnt ---------------	num of ids
		map --------------- [ids 1:cnt']


Function 5.27 neighbor_pixels: find the neighbor locations for a given foreground
Input: 	fg_map ------------ the map indicating the area or a vector indicating the locations
		im_sz ------------- size of the image [h,w] or [h,w,z]
		gap --------------- how many circles we want, 1: the most close neighbors, 2: two circles       					neighbors... default is eight neighbors
Output: neighbor_cells ---- cells contain the location of each circle

