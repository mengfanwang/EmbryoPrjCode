This file contains basic instructions for running the Matlab scripts created to import cell lineaging from the XML output files generated by the TGMM segmentation and tracking software.

1.-CONTENTS OF THE SOFTWARE ARCHIVE
-----------------------------------------------------------------------------------
We assume that the user has cloned all the files in a folder of their choice, referred to here as $ROOT_TGMM_MATLAB. This root folder contains the license agreement, this README file and multiple MATLAB scripts. 


2.-COMPILE MEX FILES
-----------------------------------------------------------------------------------

Some Matlab scripts use mex files to speed reading XML files. Execute the scipt compileMex.m to generate the required mex files.

The mex files have been tested in Linux, MacOS and Windows systems to guarantee cross-platform performance.


3.-IMPORT TRACKING RESULTS TO MATLAB ARRAY
-----------------------------------------------------------------------------------

Use the Matlab function [trackingMatrix, svIdxCell]= parseMixtureGaussiansXml2trackingMatrixCATMAIDformat(basename,frameIni,frameEnd);

INPUT:

basename:		string containing the path (except for the time point information) to the XML files generated by TGMM code. For exmaple, if they are located at E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht, then basename = 'E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht\GMEMfinalResult_frame'
frameIni:		integer with the initial time point to import (typically it is 0)
frameEnd:		integer with the final time point to import

OUTPUT:

svIdxCell:		cell array of length N, where N is the total number of objects tracked. svIdxCell{i} contains the indexes of the supervoxels belonging to the i-th object in trackingMatrix array. This index is necessary to rescue the segmentation from the .svb files output by TGMM software. The index starts in 0 following C convention. 
trackingMatrix:		numericall array of size Nx10, where N is the number of points tracked bt TGMM over time. Each of the columns contains the following information:

1.	Unique Id from the database to identify the point ( a large integer number)
2.	Cell type (represented by an integer). It is 0 if no cell type has been set for this object.
3.	x location of the nucleus centroid in world coordinates. Use the variable stackRes to convert from world coordinates to pixel unites.
4.	Same as 3 but for y location.
5.	Same as 3 but for z location.
6.	Estimated radius of the nucleus. It is 0 if this parameter was not estimated.
7.	Id of the cell in the previous time point. It is -1 if there is no linkage. Otherwise it has the unique id of the parent from column 1, so you can reconstruct the lineage.
8.	Time point of the nucleus.
9.	Confidence level in the tracking result. Value of 3 indicates high confidence that the object was correctly tracked. Value of 0 indicates low confidence.
10.	Skeleton id. All cells belonging to the same lineage have the same unique skeleton id.


4.-IMPORT TRACKING RESULTS FOR A SINGLE TIME POINT TO MATLAB STRUCTURE
-----------------------------------------------------------------------------------

Use the Matlab function obj = readXMLmixtureGaussians(filenameXML);;

INPUT:

filenameXML:		string containing the path to the XML file generated by TGMM code. For exmaple, if the TGMM oputput is located at E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht and you want to read time point 100, then filenameXML = 'E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht\GMEMfinalResult_frame0100.xml'

OUTPUT:

obj:			Matlab structure of length N, where N is the number of objects detected in the time point specified by filenameXML. For each object, the structure obj(i) contains the following fields:

-id [integer]: 			unique id of the object in this particular time point. Following C convention, the index of the first element is 0.
-lineage [integer]: 		unique id of the cell lineage the object belongs to.
-parent [integer]: 		id of the linked object at the previous time point. Following the chain of �parent� objects reconstructs the track. A value of -1 indicates the birth of a track. Following C convention, the index of the first element is 0.
-splitScore [float]: 		confidence level for the correct tracking of this particular object. A value of 0 indicates very low confidence and a value of 5 indicates very high confidence. Sorting elements by confidence level can guide the user in the data curation process and facilitate more effective editing of the TGMM results (see main text and Fig. 4).
-scale [double[3]]: 		voxel scaling factors along the x-, y- and z-axis.
-nu, beta, alpha [float]:	value of the hyper-parameters for the Bayesian GMM.
-m [double[3]]: 		mean of the Gaussian Mixture (object centroid, in pixels).
-W [double[3][3]]: 		precision matrix of the Gaussian Mixture (object shape).
-*Prior: 			same as before, but for prior values obtained from the previous time point. These values are used during the inference procedure.
-svIdx [integer[]]: 		list of indices of the super-voxels clustered by this Gaussian. Together with the �.svb� file, this information can be used to obtain precise segmentation regions for each object. Following C convention, the index of the first element is 0.
-dims [integer]:		number of dimensions of the image in each time point. This value should be 2 or 3.

5.-IMPORT SEGMENTATION INFORMATION FROM .svb FILES
-----------------------------------------------------------------------------------

Each XML file output by TGMM has a file with the same name but with the .svb extenstion (from Super-Voxel Binary file). This a custom-created binary file which stores all the indexes in the image for each super-voxel at a specific time point. In this way, and using svIdx we can recover the precise segmentation of each tracked object. 
In order to do that, you need to to use the Matlab function [svList, sizeIm] = readListSupervoxelsFromBinaryFile( filename )

INPUT:

filename:		string containing the path to the svb file generated by TGMM code. For exmaple, if the TGMM oputput is located at E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht and you want to read time point 100, then filenameXML = 'E:\TGMMruns\GMEMtracking3D_2014_6_3_15_39_21\XML_finalResult_lht\GMEMfinalResult_frame0100.svb'

OUTPUT:

svList:		cell array of length S, where S is the number of super-voxels generated at the time point specified by filename. svList{i} contains the indexes of all the voxels in the image belonging to the i-th supervoxel in ascending order. This list is exactly the same as the variable PixelIdxList returned by Matlab functions such as regionprops or bwconncomp.
		Use svList{svIdx+1} in order to access the correct super-voxel, since svIdx follow C-indexing convention.

sizeIm:		1x3 integer vector containing the size of the image stack along each dimension. This is useful to recover (x,y,z) positions from svList indexes using the function ind2sub.