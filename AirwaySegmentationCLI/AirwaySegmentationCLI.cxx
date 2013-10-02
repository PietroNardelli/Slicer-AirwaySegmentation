/*========================================================================

  Program:   Slicer
  Language:  C++
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/

#include "itkMaskNegatedImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPluginFilterWatcher.h"
#include "itkMergeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include <vector>

#include "AirwaySegmentationCLICLP.h"

#define DIM 3	

typedef signed short    InputPixelType;
typedef unsigned short  OutputPixelType;

typedef itk::Image<InputPixelType, DIM>  InputImageType;
typedef itk::Image<OutputPixelType, DIM> OutputImageType;


/** FUNCTION FOR TRACHEA SEGMENTATION */

OutputImageType::Pointer TracheaSegmentation( InputImageType::Pointer VOI, InputImageType::IndexType indexFiducialSlice, std::vector<std::vector<float> > fiducial, int labelColor = 2 )
{
	OutputImageType::Pointer trachea 	= OutputImageType::New(); 
	OutputImageType::Pointer tracheaPrev 	= OutputImageType::New();
	
	trachea->SetRegions( VOI->GetRequestedRegion() );                                   
  	trachea->SetBufferedRegion( VOI->GetBufferedRegion() );
  	trachea->SetLargestPossibleRegion( VOI->GetLargestPossibleRegion() );
  	trachea->CopyInformation( VOI );
  	trachea->Allocate();
  		
	/** TRACHEA SEGMENTATION PIPELINE */

        typedef itk::ConnectedThresholdImageFilter< InputImageType, InputImageType > ConnectedFilterType; 
	ConnectedFilterType::Pointer thresholdConnected = ConnectedFilterType::New(); 
	
	thresholdConnected->SetInput( VOI );	
	thresholdConnected->SetReplaceValue( labelColor );    
   	
	// Starting upper threshold value
	InputPixelType UpperThreshold = -850;	                    		   

	thresholdConnected->SetUpper( UpperThreshold );          
	
	InputImageType::PointType lpsPoint;
	InputImageType::IndexType index;
	
	// Seeds come in ras, convert to lps
	for( ::size_t i = 0; i < fiducial.size(); ++i )
      	{
     		lpsPoint[0] = -fiducial[i][0];
      		lpsPoint[1] = -fiducial[i][1];
     		lpsPoint[2] = fiducial[i][2];

      		VOI->TransformPhysicalPointToIndex(lpsPoint, index);
     		thresholdConnected->AddSeed( index );
	}
	         
	typedef itk::CastImageFilter<InputImageType, OutputImageType> CastingFilterType;               
       	CastingFilterType::Pointer  caster = CastingFilterType::New();		
	
	caster->SetInput( thresholdConnected->GetOutput() );  
	caster->Update();
	trachea = caster->GetOutput();	

	/** COMPUTING THE LABEL SIZES */			 
	                  
	// Extracting the axial slice containing the trachea fiducial point
 	OutputImageType::SizeType oneSliceSize;
  	
	oneSliceSize[0] = trachea->GetLargestPossibleRegion().GetSize(0);
  	oneSliceSize[1] = trachea->GetLargestPossibleRegion().GetSize(1);
  	oneSliceSize[2] = 1;
	
  	OutputImageType::RegionType axialSlice;
  	
	axialSlice.SetSize( oneSliceSize );
  	axialSlice.SetIndex( indexFiducialSlice );
 
	typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > ROIFilterType;
  	ROIFilterType::Pointer axialTracheaFilter = ROIFilterType::New();
  	
	axialTracheaFilter->SetInput( trachea );
  	axialTracheaFilter->SetRegionOfInterest( axialSlice );
  	axialTracheaFilter->Update();   

	typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > ShapeLabelType;	
  	ShapeLabelType::Pointer labelSizeFilter = ShapeLabelType::New();

  	labelSizeFilter->SetInputForegroundValue( labelColor );
  	labelSizeFilter->SetFullyConnected(1);
  	labelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
  	labelSizeFilter->Update();
	
	// Extracting the coronal slice containing the trachea fiducial point
	ROIFilterType::Pointer coronalTracheaFilter = ROIFilterType::New();
	
  	oneSliceSize[1] = 1;
  	oneSliceSize[2] = 6;
  
  	indexFiducialSlice[1] = index[1];
  	indexFiducialSlice[2] = indexFiducialSlice[2] - 3;
  
 	OutputImageType::RegionType coronalSlice;
  	coronalSlice.SetIndex( indexFiducialSlice );
 	coronalSlice.SetSize( oneSliceSize );

  	coronalTracheaFilter->SetInput( trachea );
  	coronalTracheaFilter->SetRegionOfInterest( coronalSlice );
 	coronalTracheaFilter->Update();
  
  	// Computing the sizes 
 	double xSize = 0;
 	double ySize = 0;
  	bool firstCheck = 0;
  	bool check = 0;
  	bool decrease = 0;

	do{
		if( labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() > 0 ){
	  		bool labelOverSize = 0; 
	  		for( unsigned int i = 0; i < labelSizeFilter->GetOutput()->GetNumberOfLabelObjects(); i++ ){
	      			if( labelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(1) > ySize ){
	         			ySize = labelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(1);
					if( ySize > trachea->GetLargestPossibleRegion().GetSize(1) * 0.35 ){
		    				UpperThreshold = UpperThreshold - 20;
		    				thresholdConnected->SetUpper( UpperThreshold ); 
						caster->SetInput( thresholdConnected->GetOutput() );  
						trachea = caster->GetOutput();			
						caster->Update();
						axialTracheaFilter->SetInput( trachea );
		    				axialTracheaFilter->SetRegionOfInterest( axialSlice );
		    				axialTracheaFilter->Update();
		    				i = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() - 1 ;    
		    				labelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
		    				labelSizeFilter->Update();
		    				decrease = 1;
		    				xSize = 0;
		    				ySize = 0;
	           				labelOverSize = 1;
					}
	   			}
       			}

       			if( labelOverSize == 0 ){
	  			labelSizeFilter->SetInput( coronalTracheaFilter->GetOutput() );
	  			labelSizeFilter->Update();	
	  			for( unsigned int i = 0; i < labelSizeFilter->GetOutput()->GetNumberOfLabelObjects(); i++ ){
	      				if( labelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(0) > xSize ){	
	         				xSize = labelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(0);
						if( xSize > trachea->GetLargestPossibleRegion().GetSize(0) * 2 / 3 ){
		    					UpperThreshold = UpperThreshold - 20;
		    					thresholdConnected->SetUpper( UpperThreshold ); 
							caster->SetInput( thresholdConnected->GetOutput() );  
							trachea = caster->GetOutput();			
							caster->Update();
		    					coronalTracheaFilter->SetInput( trachea );
		    					coronalTracheaFilter->SetRegionOfInterest( coronalSlice );
		    					coronalTracheaFilter->Update();
		    					i = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() - 1;	 
		    					labelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
		    					labelSizeFilter->Update();
		    					decrease = 1;
		    					xSize = 0;
		    					ySize = 0;
	        					}
             				}
         			}   
       			}
       			if( xSize != 0 && ySize != 0 ){
	  			xSize = xSize + xSize * 30 / 100;
	  			ySize = ySize + ySize * 30 / 100;
	  			firstCheck = 1;
	  		}
       		}
		else{
	  		UpperThreshold = UpperThreshold + 50;
	  		thresholdConnected->SetUpper( UpperThreshold ); 
			caster->SetInput( thresholdConnected->GetOutput() );  
			trachea = caster->GetOutput();			
			caster->Update();
	  		axialTracheaFilter->SetInput( trachea );
	  		axialTracheaFilter->SetRegionOfInterest( axialSlice );	 
	  		axialTracheaFilter->Update();
	  		coronalTracheaFilter->SetInput( trachea );
	  		coronalTracheaFilter->SetRegionOfInterest( coronalSlice );
	  		coronalTracheaFilter->Update();
		 
	  		labelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
	  		labelSizeFilter->Update();
	  		xSize = 0;
	  		ySize = 0;		
      		}
  	}
	while( !firstCheck );

	tracheaPrev = trachea;

  	/** INCREASING THE THRESHOLD ITERATIVELY UNTIL LEAKAGE OCCURS */

	typedef itk::SubtractImageFilter< OutputImageType,OutputImageType,OutputImageType > SubtractLabelImageType; 
  	SubtractLabelImageType::Pointer addedLabel = SubtractLabelImageType::New();
  
  	bool overThreshold = 0;
	
  	do{
       		addedLabel->SetInput1( trachea );
       		addedLabel->SetInput2( tracheaPrev );
       		addedLabel->Update();
       		labelSizeFilter->SetInput( addedLabel->GetOutput() );
       		labelSizeFilter->SetInputForegroundValue( labelColor );
       		labelSizeFilter->Update();
       		if( labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() > 0 ){
          		for( unsigned int i = 0; i < labelSizeFilter->GetOutput()->GetNumberOfLabelObjects(); i++ ){
	      			if( labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(0) > xSize || labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(1) > ySize || labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(2) > trachea->GetLargestPossibleRegion().GetSize(2) / 3 ){
		 			check = 1;
                 			i= labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() - 1;
	      			}
          		}
       		}
       		if( check == 0 ){
          		if( UpperThreshold < -700 ){
             			tracheaPrev = trachea;
             			if( decrease == 0 ){
                			UpperThreshold = UpperThreshold + 50;
             			}
	     			else{
	        			UpperThreshold = UpperThreshold + 10;
	     			}
	     			thresholdConnected->SetUpper( UpperThreshold ); 
				caster->SetInput( thresholdConnected->GetOutput() );  
				trachea = caster->GetOutput();			
				caster->Update();
	  		}
	  		else{
	     			check = 1;
	     			overThreshold = 1;
	  		}
      		}
  	}
  	while( check == 0 );
  
  	// Decreasing the threshold to find a better segmentation

  	if( overThreshold == 0 && decrease == 0 ){
     		while( check != 0 ){
	   		UpperThreshold = UpperThreshold - 10;
	   		thresholdConnected->SetUpper( UpperThreshold ); 
			caster->SetInput( thresholdConnected->GetOutput() );  
			trachea = caster->GetOutput();			
			caster->Update();
	   		addedLabel->SetInput1( trachea );
	   		addedLabel->SetInput2( tracheaPrev );
	   		addedLabel->Update();
	   		labelSizeFilter->SetInput( addedLabel->GetOutput() );
	   		labelSizeFilter->Update();
	   		unsigned int count = 0;
	   		if( labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() > 0 ){
	      			for( unsigned int i = 0; i < labelSizeFilter->GetOutput()->GetNumberOfLabelObjects(); i++ ){
		  			if( labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(0) < xSize && labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(1) < ySize && labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(2) < trachea->GetLargestPossibleRegion().GetSize(2) / 3){
		     				count++;
	          			}
	      			}
	      			if( count == labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() ){
		 			check = 0;
	       			}
	    		}
	    		else{
	      			check = 0;
	    		}	
     		}
  	}

	return trachea;		
}


/** FUNCTION FOR RIGHT AND LEFT AIRWAYS SEGMENTATION */

OutputImageType::Pointer RightLeftSegmentation( InputImageType::Pointer VOI, InputImageType::IndexType index, int trachea_voxels, int labelColor = 2 )
{

/** The method for the segmentation of right and left airways is based on Tschirren (2009):
	
    Tschirren, J. et al. "Airway segmentation framework for clinical environments." Proc. of Second International Workshop on Pulmonary Image Analysis. 2009.

As an alternative, Gao's method (2011) may also be used:

    Gao, D. et al. "MGRG-morphological gradient based 3D region growing algorithm for airway tree segmentation in image guided intervention therapy." Bioelectronics and Bioinformatics (ISBB), 2011 International Symposium on. IEEE, 2011.

*/
	int n_voxels 			= 0;
        int n_voxels_prev 		= 0;
        int n_voxels_max 		= trachea_voxels*0.17;
        double g;
        double g_max 			= 1.6; // 0.15 according to Gao's idea
	
	if( n_voxels_max < 15000 )
	{
		n_voxels_max = 15000;
	}


	/** SEGMENTATION PIPELINE */


	typedef itk::ConnectedThresholdImageFilter< InputImageType, InputImageType > ConnectedFilterType; 
	ConnectedFilterType::Pointer thresholdConnected = ConnectedFilterType::New();

	thresholdConnected->SetInput( VOI );			  
	thresholdConnected->SetReplaceValue( labelColor ); 
	InputPixelType UpperThreshold 	= -950;	
	thresholdConnected->SetUpper( UpperThreshold ); 
	thresholdConnected->AddSeed( index );

	typedef itk::CastImageFilter<InputImageType, OutputImageType> CastingFilterType;
	CastingFilterType::Pointer caster = CastingFilterType::New();	

	caster->SetInput( thresholdConnected->GetOutput() );  
	caster->Update();      

	// The number of voxels resulting from the first segmentation are counted

	typedef itk::StatisticsImageFilter< OutputImageType > StatisticsImageFilterType;
	StatisticsImageFilterType::Pointer StatisticsFilter = StatisticsImageFilterType::New();
	
	StatisticsFilter->SetInput(caster->GetOutput());
	StatisticsFilter->Update();
	n_voxels = StatisticsFilter->GetSum();
	

	// If the starting threshold gives an empty segmentation, neighbours are checked 

	InputImageType::SizeType radius,regionSize;
	InputImageType::IndexType regionIndex;
	InputImageType::RegionType region;	                                                        
	
		
	if( n_voxels == 0 ){
		
		bool isMinor = 0;
		
		regionSize.Fill(1);
  		regionIndex = index;
		radius.Fill(1);

  		region.SetSize(regionSize);
  		region.SetIndex(regionIndex);
		
		typedef itk::ConstNeighborhoodIterator< InputImageType > NeighborhoodIterator;
		NeighborhoodIterator iterator(radius, VOI, region);
		
		unsigned int counter = 0;

		while( counter < iterator.Size() && !isMinor )
      		{
			if( iterator.GetPixel(counter) < UpperThreshold )
			{
				index = iterator.GetIndex( counter );				
				
				thresholdConnected->ClearSeeds();
				thresholdConnected->AddSeed( index );
				thresholdConnected->Update();
				
				caster->SetInput( thresholdConnected->GetOutput() );   	
				caster->Update();	
			
				isMinor = 1;
 			}
			counter++;   		
    		}
		if ( !isMinor )
		{
	  		std::cout<<"Please move the seed point in a different position."<<std::endl;
			return caster->GetOutput();
		}

	StatisticsFilter->SetInput(caster->GetOutput());
	StatisticsFilter->Update();
	n_voxels = StatisticsFilter->GetSum();
	}
	
	// If the number of voxels resulting form the segmentation is too high the threshold is iteratively decreased

	if( n_voxels > n_voxels_max ){
		
		while( n_voxels > n_voxels_max ){
			
			UpperThreshold = UpperThreshold - 20;
			thresholdConnected->SetUpper( UpperThreshold );
			
			caster->SetInput( thresholdConnected->GetOutput() );  
			caster->Update();
				
			StatisticsFilter->SetInput(caster->GetOutput());
			StatisticsFilter->Update();			
			n_voxels = StatisticsFilter->GetSum();
				
			if( n_voxels == 0 ){
				bool isMinor = 0;
		
				regionIndex = index;

  				regionSize.Fill(1);	
				radius.Fill(1);

				region.SetSize(regionSize);
  				region.SetIndex(regionIndex);
				
				typedef itk::ConstNeighborhoodIterator< InputImageType > NeighborhoodIterator;
				NeighborhoodIterator iterator(radius, VOI, region);
				
				unsigned int counter = 0;
				
				while( counter < iterator.Size() && !isMinor)
      				{
					if( iterator.GetPixel(counter) < UpperThreshold )
					{
						index = iterator.GetIndex( counter );				
						
						thresholdConnected->ClearSeeds();
						thresholdConnected->AddSeed( index );
						thresholdConnected->Update(); 

						caster->SetInput( thresholdConnected->GetOutput() );   	
						caster->Update();
				
						isMinor = 1;		
 					}
					counter++;   		
    				}
				StatisticsFilter->SetInput(caster->GetOutput());
				StatisticsFilter->Update();
				n_voxels = StatisticsFilter->GetSum();		
			}
		}

		// If n_voxels is too small, an increase of the threshold of even 1 HU might cause the violation of g < g_max

		while( n_voxels < 3000 ){ 
			UpperThreshold = UpperThreshold + 1;
			thresholdConnected->SetUpper( UpperThreshold );

			caster->SetInput( thresholdConnected->GetOutput() );  
			caster->Update();
			
			StatisticsFilter->SetInput(caster->GetOutput());
			StatisticsFilter->Update();	
			
			n_voxels = StatisticsFilter->GetSum();
		}
			
		do{	
			UpperThreshold = UpperThreshold + 1;
			n_voxels_prev = n_voxels;				
			
			thresholdConnected->SetUpper( UpperThreshold );	
			
			caster->SetInput( thresholdConnected->GetOutput() );  
			caster->Update();
	
			StatisticsFilter->SetInput(caster->GetOutput());
			StatisticsFilter->Update();	

			n_voxels = StatisticsFilter->GetSum();
			g = double( n_voxels/n_voxels_prev );	// double((n_voxels - n_voxels_prev)/n_voxels_prev) according to Gao et al.
		}while( g < g_max && n_voxels < n_voxels_max );
			
		UpperThreshold = UpperThreshold - 1;
		thresholdConnected->SetUpper( UpperThreshold );

		caster->SetInput( thresholdConnected->GetOutput() ); 
		caster->Update();
	}

	// The threshold is iteratively increased until leakage occurs

	else{

		// If n_voxels is too small, an increase of the threshold of even 1 HU might cause the violation of g < g_max		
		
		while( n_voxels < 3000 ){ 
			UpperThreshold = UpperThreshold + 1;
			thresholdConnected->SetUpper( UpperThreshold );

			caster->SetInput( thresholdConnected->GetOutput() );  
			caster->Update();
			
			StatisticsFilter->SetInput(caster->GetOutput());
			StatisticsFilter->Update();	
			
			n_voxels = StatisticsFilter->GetSum();
		}
		do{	
			UpperThreshold = UpperThreshold + 20;
			n_voxels_prev = n_voxels;	
			
			thresholdConnected->SetUpper( UpperThreshold );
	
			caster->SetInput( thresholdConnected->GetOutput() );  
			caster->Update();
			
			StatisticsFilter->SetInput(caster->GetOutput());
			StatisticsFilter->Update();

			n_voxels = StatisticsFilter->GetSum();
			g = double( n_voxels/n_voxels_prev );	// double((n_voxels - n_voxels_prev)/n_voxels_prev) according to Gao et al.
		}while( g < g_max && n_voxels < n_voxels_max );
		
		UpperThreshold = UpperThreshold - 20;		
		thresholdConnected->SetUpper( UpperThreshold );

		caster->SetInput( thresholdConnected->GetOutput() );  
		caster->Update();
	
		StatisticsFilter->SetInput(caster->GetOutput());
		StatisticsFilter->Update();	
		n_voxels = StatisticsFilter->GetSum();
		
		do{	
			UpperThreshold = UpperThreshold + 1;
			n_voxels_prev = n_voxels;	
			
			thresholdConnected->SetUpper( UpperThreshold );		
	
			caster->SetInput( thresholdConnected->GetOutput() );  
			caster->Update();

			StatisticsFilter->SetInput(caster->GetOutput());
			StatisticsFilter->Update();
			
			n_voxels = StatisticsFilter->GetSum();
			g = double( n_voxels/n_voxels_prev );	// double((n_voxels - n_voxels_prev)/n_voxels_prev) according to Gao et al.
		}while( g < g_max && n_voxels < n_voxels_max );
		
		UpperThreshold = UpperThreshold - 1;
		thresholdConnected->SetUpper( UpperThreshold );	
		
		caster->SetInput( thresholdConnected->GetOutput() );  
		caster->Update();
	}
	
	StatisticsFilter->SetInput(caster->GetOutput());
	StatisticsFilter->Update();			
	n_voxels = StatisticsFilter->GetSum();
	
	return caster->GetOutput(); 
}


/** FUNCTION WHICH PASTES AN IMAGE IN A SPECIFIED INDEX OF THE DESTINATION IMAGE */

template <class ImageType>													                     
typename ImageType::Pointer Paste( typename ImageType::Pointer sourceImage, typename ImageType::IndexType index, typename ImageType::Pointer destImage)
{
	typedef itk::PasteImageFilter< ImageType, ImageType > PasteImageFilterType;
	typename PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New ();  
      
	pasteFilter->SetSourceImage(sourceImage);                                                 
	pasteFilter->SetDestinationImage(destImage);                                               
	pasteFilter->SetSourceRegion(sourceImage->GetLargestPossibleRegion());                    
	pasteFilter->SetDestinationIndex(index);                                                  

	try
	{
		pasteFilter->Update();				                 	   
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	} 
		
	return pasteFilter->GetOutput();
}



int main( int argc, char *argv[] )
{

	PARSE_ARGS; 	                         
  	
	typedef  itk::ImageFileReader<InputImageType>  ReaderType;				
  	ReaderType::Pointer reader = ReaderType::New();     
  	reader->SetFileName(  inputVolume.c_str() );                                            
  	try
    	{
		reader->Update();
    	}
	catch ( itk::ExceptionObject & e )
    	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		return EXIT_FAILURE;
    	}

  	// The different labels will be pasted into pasteImage

	OutputImageType::Pointer pasteImage = OutputImageType::New(); 

  	pasteImage->SetRegions( reader->GetOutput()->GetRequestedRegion() );                                   
  	pasteImage->SetBufferedRegion( reader->GetOutput()->GetBufferedRegion() );
  	pasteImage->SetLargestPossibleRegion( reader->GetOutput()->GetLargestPossibleRegion() );
  	pasteImage->CopyInformation( reader->GetOutput() );
  	pasteImage->Allocate();
  	pasteImage->FillBuffer(0);
	
  	// Seed points indices

  	InputImageType::IndexType tracheaFiducial, rightFiducial, leftFiducial;
  	InputImageType::PointType tracheaPoint, rightPoint, leftPoint;

  	std::vector< std::vector<float> > tracheaSeedPoint, rightSeedsVector, leftSeedsVector;
  
	// Finding the trachea seed point 

	int tracheaIndex = 0;

  	if( seed.size() > 0 && seed.size() < 4 ){ 					
   	 	
   	 	for( ::size_t i = 0; i < seed.size(); ++i ){
     	 		if( seed[i][2] > seed[tracheaIndex][2] ){
	   			tracheaIndex = i;
			}
    		}       
    		tracheaSeedPoint.push_back( seed[tracheaIndex] );

    		// Convert to lps the seed point
    		tracheaPoint[0] = -seed[tracheaIndex][0];	
   		tracheaPoint[1] = -seed[tracheaIndex][1];
    		tracheaPoint[2] =  seed[tracheaIndex][2];

    		// Convert the lps physical point to index
    		reader->GetOutput()->TransformPhysicalPointToIndex( tracheaPoint, tracheaFiducial );	     		
  	}
  	else{
    		if( seed.size() == 0 ){
			std::cerr << "No seeds specified!" << std::endl;
    			return -1;
		}
		else{
			std::cerr << "Please place only three seed points!" << std::endl;
			return -1;
		} 
  	}   
  
  	/** TRACHEA SEGMENTATION */

  	InputImageType::SizeType cropSize;	                    
  	InputImageType::IndexType tracheaCropStart;              
  
 	if( seed.size() > 1 ){

		// Main bronchi seed points	
   		int prevRight, prevLeft;
 		prevRight = tracheaIndex;
    		prevLeft  = tracheaIndex; 
    		
		int rightPosition, leftPosition; 
	
    		for( ::size_t i = 0; i < seed.size(); ++i ){
			if( seed[i][0] > seed[tracheaIndex][0] ){         // Is the seed on the right side of the trachea point...	
	  			rightSeedsVector.push_back( seed[i] );
	  			if( seed[i][0] < seed[prevRight][0] ){
	     				rightPosition = i;
	     				prevRight = rightPosition;
	  			}
	  			else if( prevRight != tracheaIndex ){
	     				rightPosition = prevRight;
	  			}
	  			else{
	     				rightPosition = i;
	      				prevRight = i;
	  			}
        		}
	  
        		if( seed[i][0] < seed[tracheaIndex][0] ){         //...or on the left one?
	  			leftSeedsVector.push_back( seed[i] );
	  			if( seed[i][0] > seed[prevLeft][0] ){
	   				leftPosition = i;
	   				prevLeft = i;
	  			}
	  			else if( prevLeft != tracheaIndex ){
	   				leftPosition = prevLeft;
	  			}
	  			else{
	   				leftPosition = i;
	   				prevLeft = i;
	  			}
        		}     
    		}
   
		// Converting the right seed point to lps
    		if( !rightSeedsVector.empty() ){			
      			rightPoint[0] = -seed[rightPosition][0];	
      			rightPoint[1] = -seed[rightPosition][1];
      			rightPoint[2] =  seed[rightPosition][2];
			reader->GetOutput()->TransformPhysicalPointToIndex( rightPoint, rightFiducial );
    		}

		// Converting the left seed point to lps
    		if( !leftSeedsVector.empty() ){
      			leftPoint[0] = -seed[leftPosition][0];		
      			leftPoint[1] = -seed[leftPosition][1];
      			leftPoint[2] =  seed[leftPosition][2];
			reader->GetOutput()->TransformPhysicalPointToIndex( leftPoint, leftFiducial );
    		}

    		if( rightSeedsVector.empty() ){  // No right seed point
      			cropSize[0] = ( leftFiducial[0] - tracheaFiducial[0] ) * 2;
      			cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - leftFiducial[2] - 1;	
			tracheaCropStart[0] = leftFiducial[0] - cropSize[0];       			
    		}
    		else if( leftSeedsVector.empty() ){  // No left seed point
      			cropSize[0] = ( tracheaFiducial[0] - rightFiducial[0] ) * 2;
      			cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - rightFiducial[2] - 1;
			tracheaCropStart[0] = rightFiducial[0] + 1;   
    		}
    		else{								 
       			cropSize[0] =  leftFiducial[0] - rightFiducial[0] - 2;
       			if( rightFiducial[2] > leftFiducial[2] ){
         			cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - rightFiducial[2] - 1;
       			}
       			else{
         			cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - leftFiducial[2] - 1;
       			}
			tracheaCropStart[0] = rightFiducial[0] + 1; 
      		}

		tracheaCropStart[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - cropSize[2] - 1; 
  	}

  	else{ // Only one fiducial was set
    		cropSize[0] = 60; // To be fixed
    		cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2);
		tracheaCropStart[0] = tracheaFiducial[0] - 30;   // To be fixed
    		tracheaCropStart[2] = 0;
  	}
 
  	cropSize[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(1);        
  	tracheaCropStart[1] = 0;		     
  
  	InputImageType::RegionType DesiredRegion;                                                  
  	DesiredRegion.SetSize(  cropSize  );                                                                
  	DesiredRegion.SetIndex( tracheaCropStart );                                                          
  
	// Cropping the trachea 
	typedef itk::RegionOfInterestImageFilter< InputImageType, InputImageType > inputROIFilterType;  
  	inputROIFilterType::Pointer ROIFilter = inputROIFilterType::New();	                 
  
  	ROIFilter->SetInput( reader->GetOutput() );						          
  	ROIFilter->SetRegionOfInterest( DesiredRegion );					 
  	ROIFilter->Update();				
    
  	/** SEGMENTING THE TRACHEA */ 

  	OutputImageType::Pointer trachea = OutputImageType::New();
	
  	InputImageType::IndexType FiducialSlice;
  	FiducialSlice.Fill(0);
  	FiducialSlice[2] = cropSize[2] - ( reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - tracheaFiducial[2] ); 

	trachea = TracheaSegmentation( ROIFilter->GetOutput(), FiducialSlice, tracheaSeedPoint, labelValue );	
    
  	OutputImageType::IndexType regionIndex;

  	regionIndex = tracheaCropStart;
  	pasteImage = Paste<OutputImageType>( trachea, regionIndex, pasteImage );   

	/** RIGHT AND LEFT AIRWAYS SEGMENTATION*/

  	if( seed.size() > 1 ){ 

		// The trachea label is used to mask the input image
		typedef itk::MaskNegatedImageFilter< InputImageType, OutputImageType, InputImageType > 	MaskNegatedImageType;
		MaskNegatedImageType::Pointer maskFilter = MaskNegatedImageType::New();     

		maskFilter->SetInput1( reader->GetOutput() );
     		maskFilter->SetInput2( pasteImage );
     		maskFilter->Update();	

     		typedef itk::StatisticsImageFilter< OutputImageType > StatisticsImageFilterType;		
 		StatisticsImageFilterType::Pointer StatisticsFilter = StatisticsImageFilterType::New();	

		StatisticsFilter->SetInput(trachea);
		StatisticsFilter->Update();
		
		typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > ShapeLabelType;
		
		typedef ShapeLabelType::OutputImageType LabelMapType;
		LabelMapType::Pointer labelMap = LabelMapType::New();
		
		/** RIGHT AIRWAY SEGMENTATION*/

     		if( !rightSeedsVector.empty() ){
			
			OutputImageType::Pointer rightLung = OutputImageType::New();

			rightLung = RightLeftSegmentation( maskFilter->GetOutput(), rightFiducial, StatisticsFilter->GetSum(), labelValue ); 
					
			ShapeLabelType::Pointer rightLabelConverter = ShapeLabelType::New();
		  
			rightLabelConverter->SetInput( rightLung );
			rightLabelConverter->SetInputForegroundValue( labelValue );
			rightLabelConverter->Update();
		
			labelMap = rightLabelConverter->GetOutput();		
     		}

		/** LEFT AIRWAY SEGMENTATION*/

     		if( !leftSeedsVector.empty() ){

			OutputImageType::Pointer leftLung = OutputImageType::New();
		
			leftLung = RightLeftSegmentation( maskFilter->GetOutput(), leftFiducial, StatisticsFilter->GetSum(), labelValue ); 
				
			ShapeLabelType::Pointer leftLabelConverter = ShapeLabelType::New();

			leftLabelConverter->SetInput( leftLung );
			leftLabelConverter->SetInputForegroundValue( labelValue );
			leftLabelConverter->Update();  
			
			if( !rightSeedsVector.empty() ){

				typedef itk::MergeLabelMapFilter< LabelMapType > MergeFilterType;	  
				MergeFilterType::Pointer mergeFilter = MergeFilterType::New(); 
		  
				mergeFilter->SetMethod( MergeFilterType::PACK );
				mergeFilter->SetInput( labelMap );
				mergeFilter->SetInput( 1, leftLabelConverter->GetOutput() );
				mergeFilter->Update();

	        		labelMap = mergeFilter->GetOutput();
			}
			else{
				labelMap = leftLabelConverter->GetOutput();
			}
     		}	

		ShapeLabelType::Pointer labelConverter = ShapeLabelType::New();

     		labelConverter->SetInput( pasteImage );
     		labelConverter->SetInputForegroundValue( labelValue ); 
     		labelConverter->Update();
	  
     		for( unsigned int i = 0; i < labelConverter->GetOutput()->GetNumberOfLabelObjects(); i++ ){
	 		labelMap->PushLabelObject( labelConverter->GetOutput()->GetNthLabelObject(i) );
     		}
	  
     		labelMap->Update();
	  
		typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > LabelMapToBinaryImageType; 
     		LabelMapToBinaryImageType::Pointer labelToBinaryFilter = LabelMapToBinaryImageType::New();

     		labelToBinaryFilter->SetInput( labelMap );
     		labelToBinaryFilter->SetBackgroundValue( 0 );
     		labelToBinaryFilter->SetForegroundValue( labelValue );
     		labelToBinaryFilter->Update();
	 
     		pasteImage = labelToBinaryFilter->GetOutput();
  	}

  	/** CLOSING AND HOLE FILLING TO IMPROVE THE SEGMENTATION */

	typedef itk::BinaryBallStructuringElement< OutputImageType::PixelType, DIM > StructuringElementType;
  	
	StructuringElementType structElement;
	StructuringElementType::SizeType radius;
  	radius.Fill( 2 );

  	structElement.SetRadius( radius );
  	structElement.CreateStructuringElement();
	
	typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > CloseType;
 	CloseType::Pointer closing = CloseType::New();

	itk::PluginFilterWatcher watcher1(closing, "Closing Operation", CLPProcessInformation);

  	closing->SetInput( pasteImage );
 	closing->SetKernel( structElement );
  	closing->SetForegroundValue( labelValue );
  	closing->Update();
	
	typedef itk::VotingBinaryIterativeHoleFillingImageFilter< OutputImageType > IterativeFillHolesFilterType;
	IterativeFillHolesFilterType::Pointer HoleFilling= IterativeFillHolesFilterType::New();

  	itk::PluginFilterWatcher watcher2(HoleFilling, "Holes Filling Operation", CLPProcessInformation);

  	OutputImageType::SizeType FillRadius;

  	FillRadius.Fill(1);

	HoleFilling->SetInput( closing->GetOutput() );
  	HoleFilling->SetRadius( FillRadius );
	HoleFilling->SetBackgroundValue( 0 );
  	HoleFilling->SetForegroundValue( labelValue );
  	HoleFilling->SetMajorityThreshold( 1 );
  	HoleFilling->SetMaximumNumberOfIterations( 15 );
	HoleFilling->Update();

  	/** LABEL CREATION */

	typedef  itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer labelImage = WriterType::New();

	labelImage->SetFileName( label.c_str() );
  	labelImage->SetInput( HoleFilling->GetOutput() );
	labelImage->SetUseCompression(1);

  	try
    	{
		labelImage->Update();
    	}
	catch ( itk::ExceptionObject & e )
    	{
		std::cerr << "exception in file writer " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		return EXIT_FAILURE;
    	}

  	return EXIT_SUCCESS;
}
