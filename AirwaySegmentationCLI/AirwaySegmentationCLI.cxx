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
#include "itkFlipImageFilter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageDuplicator.h"

#include "itkGrayscaleFillholeImageFilter.h"

#include <vector>

#include "AirwaySegmentationCLICLP.h"

#define DIM 3


typedef signed short    InputPixelType;
typedef unsigned short  OutputPixelType;

typedef itk::Image<InputPixelType, DIM>  InputImageType;
typedef itk::Image<OutputPixelType, DIM> OutputImageType;

/** FUNCTION FOR TRACHEA SEGMENTATION */

OutputImageType::Pointer TracheaSegmentation( InputImageType::Pointer VOI,
                                              InputImageType::IndexType indexFiducialSlice,
                                              std::vector<std::vector<float> > fiducial,
                                              int labelColor = 2 )
{
    OutputImageType::Pointer trachea 		= OutputImageType::New(); 
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
    InputPixelType UpperThreshold = -900;	                    		   

    thresholdConnected->SetUpper( UpperThreshold );          
	
    InputImageType::PointType lpsPoint;
    InputImageType::IndexType index;
	
    // Seeds come in ras, convert to lps
    for( ::size_t i = 0; i < fiducial.size(); ++i )
    {
        lpsPoint[0] = fiducial[i][0] * (-VOI->GetDirection()[0][0]);
        lpsPoint[1] = fiducial[i][1] * (-VOI->GetDirection()[1][1]);
        lpsPoint[2] = fiducial[i][2] *   VOI->GetDirection()[2][2];

        VOI->TransformPhysicalPointToIndex(lpsPoint, index);
        thresholdConnected->AddSeed( index );
    }

    typedef itk::CastImageFilter<InputImageType, OutputImageType> CastingFilterType;  
    CastingFilterType::Pointer  caster = CastingFilterType::New();		
	
    caster->SetInput( thresholdConnected->GetOutput() );  
    caster->Update();
    trachea = caster->GetOutput();	

    /** COMPUTING THE LABEL SIZES */	 	                  
    OutputImageType::Pointer tracheaAxialCopy = OutputImageType::New(); 
    OutputImageType::Pointer tracheaCoronalCopy = OutputImageType::New(); 

    typedef itk::ImageDuplicator<OutputImageType> DuplicatorFilterType;

    DuplicatorFilterType::Pointer duplicatorFilter = DuplicatorFilterType::New();
    duplicatorFilter->SetInputImage(trachea);
    duplicatorFilter->Update();

    // Extracting the axial slice containing the trachea fiducial point
    OutputImageType::SizeType  oneAxialSliceSize;
    InputImageType::IndexType  indexAxialSlice = indexFiducialSlice;
  	
    oneAxialSliceSize[0] = trachea->GetLargestPossibleRegion().GetSize(0);
    oneAxialSliceSize[1] = trachea->GetLargestPossibleRegion().GetSize(1);
    unsigned int diff = trachea->GetLargestPossibleRegion().GetSize(2)-indexAxialSlice[2];
    if( trachea->GetLargestPossibleRegion().GetSize(2) > 40 &&
        indexAxialSlice[2] >= 20 &&
        diff >= 20 )
    {
        oneAxialSliceSize[2] = 40;
        indexAxialSlice[2]  -= 20;
    }
    else if( trachea->GetLargestPossibleRegion().GetSize(2) > 40 &&
             indexAxialSlice[2] >= 20 &&
             diff < 20 )
    {
        oneAxialSliceSize[2] = 40;
        indexAxialSlice[2]   = trachea->GetLargestPossibleRegion().GetSize(2) - 40;
    }
    else if( trachea->GetLargestPossibleRegion().GetSize(2) > 40 && indexAxialSlice[2] < 20 )
    {
        oneAxialSliceSize[2] = 40;
        indexAxialSlice  [2] = 0;
    }
    else if( trachea->GetLargestPossibleRegion().GetSize(2) <= 40 )
    {
        oneAxialSliceSize[2] = trachea->GetLargestPossibleRegion().GetSize(2);
        indexAxialSlice  [2] = 0;
    }
	
    OutputImageType::RegionType axialSlice;
 
    typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > ROIFilterType;
    ROIFilterType::Pointer axialTracheaFilter = ROIFilterType::New();

    typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > ShapeLabelType;	
    ShapeLabelType::Pointer axialLabelSizeFilter = ShapeLabelType::New();

    axialLabelSizeFilter->SetInputForegroundValue( labelColor );
    axialLabelSizeFilter->SetFullyConnected(1);
	
    // Extracting the coronal slice containing the trachea fiducial point
    OutputImageType::SizeType oneCoronalSliceSize;
    oneCoronalSliceSize[0] = trachea->GetLargestPossibleRegion().GetSize(0);
    oneCoronalSliceSize[1] = 1;
    oneCoronalSliceSize[2] = 6;
  
    InputImageType::IndexType indexCoronalSlice;
    indexCoronalSlice.Fill(0);
    indexCoronalSlice[1] = index[1];
    if( indexFiducialSlice[2] >= 3 )
    {
        indexCoronalSlice[2] = indexFiducialSlice[2] - 3;
    }
    else
    {
        indexCoronalSlice[2] = indexFiducialSlice[2];
    }
    OutputImageType::RegionType coronalSlice;

    ROIFilterType::Pointer coronalTracheaFilter = ROIFilterType::New();

    ShapeLabelType::Pointer coronalLabelSizeFilter = ShapeLabelType::New();
    coronalLabelSizeFilter->SetInputForegroundValue( labelColor );
    coronalLabelSizeFilter->SetFullyConnected(1);
  
    // Computing the sizes 
    double xSize      = 0;
    double ySize      = 0;
    bool   firstCheck = 0;
    bool   check      = 0;
    bool   decrease   = 0;
       
    double tracheaYSize  = trachea->GetLargestPossibleRegion().GetSize(1) * 0.25;
    double tracheaXSize  = trachea->GetLargestPossibleRegion().GetSize(0) * 2/3;

    do{
        axialSlice.SetSize( oneAxialSliceSize );
        axialSlice.SetIndex( indexAxialSlice );

        duplicatorFilter->Update();	
        tracheaAxialCopy = duplicatorFilter->GetOutput();	
            
        axialTracheaFilter->SetInput( tracheaAxialCopy );
        axialTracheaFilter->SetRegionOfInterest( axialSlice );
        axialTracheaFilter->Update();

        axialLabelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
        axialLabelSizeFilter->Update();

        if( axialLabelSizeFilter->GetOutput()->GetNumberOfLabelObjects() > 0 )
        {
            bool labelOverSize = 0; 
            unsigned int numberOfObjects = axialLabelSizeFilter->GetOutput()->GetNumberOfLabelObjects();

            for( unsigned int i = 0; i < numberOfObjects; ++i )
            {
                ySize += axialLabelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(1);
                if( ySize > tracheaYSize)
                {
                    UpperThreshold = UpperThreshold - 20;

                    thresholdConnected->SetUpper( UpperThreshold ); 
                    caster->SetInput( thresholdConnected->GetOutput() ); 
                    caster->Update();
                    trachea = caster->GetOutput();
                    duplicatorFilter->Update();	
                    tracheaAxialCopy = duplicatorFilter->GetOutput();	
                    axialTracheaFilter->SetInput( tracheaAxialCopy );
                    axialTracheaFilter->SetRegionOfInterest( axialSlice );
                    axialTracheaFilter->Update();
                    axialLabelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
                    axialLabelSizeFilter->Update();
                    decrease = 1;
                    xSize = 0;
                    ySize = 0;
                    labelOverSize = 1;
                    i = numberOfObjects - 1;    
                }
            }

            if( !labelOverSize )
            {
                coronalSlice.SetIndex( indexCoronalSlice );
                coronalSlice.SetSize( oneCoronalSliceSize );

                duplicatorFilter->Update();	
                tracheaCoronalCopy = duplicatorFilter->GetOutput();	

                coronalTracheaFilter->SetInput( tracheaCoronalCopy );
                coronalTracheaFilter->SetRegionOfInterest( coronalSlice );
                coronalTracheaFilter->Update();

                coronalLabelSizeFilter->SetInput( coronalTracheaFilter->GetOutput() );
                coronalLabelSizeFilter->Update();
		 
                unsigned int numberOfObjects = coronalLabelSizeFilter->GetOutput()->GetNumberOfLabelObjects();

                for( unsigned int i = 0; i < numberOfObjects; i++ )
                {
                    xSize += coronalLabelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(0);

                    if( xSize > tracheaXSize )
                    {
                        UpperThreshold = UpperThreshold - 20;

                        thresholdConnected->SetUpper( UpperThreshold ); 
                        caster->SetInput( thresholdConnected->GetOutput() );  
                        caster->Update();
                        trachea = caster->GetOutput();
                        duplicatorFilter->Update();	
                        tracheaCoronalCopy = duplicatorFilter->GetOutput();	
                        coronalTracheaFilter->SetInput( tracheaCoronalCopy );
                        coronalTracheaFilter->SetRegionOfInterest( coronalSlice );
                        coronalTracheaFilter->Update();
                        i = numberOfObjects - 1;	 
                        coronalLabelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
                        coronalLabelSizeFilter->Update();
                        decrease = 1;
                        xSize = 0;
                        ySize = 0;
                    }
                }
            }
            if( xSize != 0 && ySize != 0 )
            {
                xSize = xSize + xSize * 30 / 100;
                ySize = ySize + ySize * 30 / 100;
                firstCheck = 1;
            }
        }
        else
        {
            bool isMinor = 0;
            InputImageType::SizeType    radius,regionSize;
            InputImageType::IndexType   regionIndex;
            InputImageType::RegionType  region;	  

            regionSize.Fill(3);
            regionIndex = index;
            radius.Fill(3);

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
            
                    indexCoronalSlice[1] = index[1];
                    indexCoronalSlice[2] = index[2] - 3;

                    thresholdConnected->ClearSeeds();
                    thresholdConnected->AddSeed( index );
                    thresholdConnected->Update();                
    
                    caster->SetInput( thresholdConnected->GetOutput() );
                    caster->Update();	

                    trachea = caster->GetOutput();	
                    isMinor = 1;
                }
                counter++;   
            }

            if ( !isMinor && !decrease)
            {
                if( UpperThreshold < -800 )
                {
                    UpperThreshold = UpperThreshold + 50;

                    thresholdConnected->SetUpper( UpperThreshold ); 
                    thresholdConnected->Update();                
    
                    caster->SetInput( thresholdConnected->GetOutput() );
                    caster->Update();

                    trachea = caster->GetOutput();
                }
                else
                {
                    std::cout<<"Please move the seed point in a different position."<<std::endl;
                    return trachea;
                }
            }
            else if( !isMinor && decrease )
            {
                if( UpperThreshold < -800 )
                {
                    UpperThreshold = UpperThreshold + 1;

                    thresholdConnected->SetUpper( UpperThreshold ); 
                    thresholdConnected->Update();

                    caster->SetInput( thresholdConnected->GetOutput() ); 
                    caster->Update();

                    trachea = caster->GetOutput();	
                }
                else
                {
                    std::cout<<"Please move the seed point to a different location."<<std::endl;
                    return trachea;
                }
            }
			
            axialSlice.SetSize( oneAxialSliceSize );
            axialSlice.SetIndex( indexAxialSlice );

            duplicatorFilter->Update();	
            tracheaAxialCopy = duplicatorFilter->GetOutput();	
            
            axialTracheaFilter->SetInput( tracheaAxialCopy );
            axialTracheaFilter->SetRegionOfInterest( axialSlice );
            axialTracheaFilter->Update();

            axialLabelSizeFilter->SetInput( axialTracheaFilter->GetOutput() );
            axialLabelSizeFilter->Update();

            coronalSlice.SetIndex( indexCoronalSlice );
            coronalSlice.SetSize( oneCoronalSliceSize );

            duplicatorFilter->Update();	
            tracheaCoronalCopy = duplicatorFilter->GetOutput();	

            coronalTracheaFilter->SetInput( tracheaCoronalCopy );
            coronalTracheaFilter->SetRegionOfInterest( coronalSlice );
            coronalTracheaFilter->Update();

            coronalLabelSizeFilter->SetInput( coronalTracheaFilter->GetOutput() );
            coronalLabelSizeFilter->Update();
	  		
            xSize = 0;
            ySize = 0;
        }
    }
    while( !firstCheck && UpperThreshold > -1100 );
    
    duplicatorFilter->SetInputImage(trachea);
    duplicatorFilter->Update();
    tracheaPrev = duplicatorFilter->GetOutput();
        
    /** INCREASING THE THRESHOLD ITERATIVELY UNTIL LEAKAGE OCCURS */
    typedef itk::SubtractImageFilter< OutputImageType,OutputImageType,OutputImageType > SubtractLabelImageType; 
    SubtractLabelImageType::Pointer addedLabel = SubtractLabelImageType::New();
  
    bool overThreshold = 0;

    ShapeLabelType::Pointer labelSizeFilter = ShapeLabelType::New();

    do{
        addedLabel->SetInput1( trachea );
        addedLabel->SetInput2( tracheaPrev );
        addedLabel->Update();

        labelSizeFilter->SetInput( addedLabel->GetOutput() );
        labelSizeFilter->SetInputForegroundValue( labelColor );
        labelSizeFilter->Update();
        unsigned int numberOfObjects = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects(); 
        double       xSz             = 0;
        double       ySz             = 0;
        double       zSz             = 0;
        
        if( numberOfObjects > 0 )
        {
            for( unsigned int i = 0; i < numberOfObjects; i++ )
            {
                xSz = labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(0);
                ySz = labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(1);
                zSz = labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(2);
                if( xSz > xSize || ySz > ySize || zSz > trachea->GetLargestPossibleRegion().GetSize(2) / 3 )
                {
                    if( decrease )
                    {
                        UpperThreshold = UpperThreshold - 10;
                        thresholdConnected->SetUpper( UpperThreshold ); 
                        thresholdConnected->Update();
                        caster->SetInput( thresholdConnected->GetOutput() );
                        caster->Update();
                        trachea = caster->GetOutput();
                    }
                    check = 1;
                    i = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() - 1;
                }
            }
        }
        if( !check )
        {
            if( UpperThreshold < -800 )
            {
                duplicatorFilter->SetInputImage(trachea);
                duplicatorFilter->Update();
                tracheaPrev = duplicatorFilter->GetOutput();
                if( !decrease )
                {
                    UpperThreshold = UpperThreshold + 50;
                }
                else
                {
                    UpperThreshold = UpperThreshold + 10;
                }
                thresholdConnected->SetUpper( UpperThreshold ); 
                thresholdConnected->Update();
                caster->SetInput( thresholdConnected->GetOutput() );
                caster->Update();
                trachea = caster->GetOutput();	
            }
            else
            {
                check = 1;
                overThreshold = 1;
            }
        }
    }
    while( !check );
  
    // Decreasing the threshold to find a better segmentation
    if( !overThreshold && !decrease )
    {
        while( check && UpperThreshold > -1100 )
        {
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
            unsigned int numberOfObjects = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects();

            if( numberOfObjects > 0 )
            {
                for( unsigned int i = 0; i < numberOfObjects; i++ )
                {
                    if( labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(0) < xSize && 
                        labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(1) < ySize &&
                        labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(2) < trachea->GetLargestPossibleRegion().GetSize(2) / 3)
                    {
                        count++;
                    }
                }
                if( count == numberOfObjects )
                {
                    check = 0;
                }
            }
            else
            {
                check = 0;
            }
        }
    }

    return trachea;	
}


/** FUNCTION FOR RIGHT AND LEFT AIRWAYS SEGMENTATION */
OutputImageType::Pointer RightLeftSegmentation( InputImageType::Pointer VOI, 
                                                InputImageType::IndexType index,
                                                std::string reconKernel,
                                                int trachea_voxels, 
                                                int labelColor = 2 )
{
    /** The method for the segmentation of right and left airways is based on Tschirren (2009):
        Tschirren, J. et al. "Airway segmentation framework for clinical environments." Proc. of Second International Workshop on Pulmonary Image Analysis. 2009.
        As an alternative, Gao's method (2011) may also be used:
        Gao, D. et al. "MGRG-morphological gradient based 3D region growing algorithm for airway tree segmentation in image guided intervention therapy." Bioelectronics and Bioinformatics (ISBB), 2011 International Symposium on. IEEE, 2011.
    */ 

    double th = 0;

    if( reconKernel == "STANDARD" || reconKernel == "B20f" || reconKernel == "B30f" || reconKernel == "B" 
        || reconKernel == "C" || reconKernel == "FC10" || reconKernel == "FC12")
    {
        if( VOI->GetLargestPossibleRegion().GetSize(2) <= 300 )
        {
            if( trachea_voxels > 50000 )
            {
                th = 0.5;
            }
            else if( trachea_voxels > 20000 && trachea_voxels <= 50000 )//( trachea_voxels > 40000 && trachea_voxels <= 100000 )
            {
                th = 0.75;
            }
            else if( trachea_voxels <= 20000 )
            {
                th = 0.9;
            }
        }
        if( VOI->GetLargestPossibleRegion().GetSize(2) > 300 && VOI->GetLargestPossibleRegion().GetSize(2) <= 400 )
        {
            if( trachea_voxels > 100000 ) 
            {
                th = 0.5;
            }
            else if( trachea_voxels > 85000 && trachea_voxels <= 100000 )
            {
                th = 0.75;
            }
            else if( trachea_voxels <= 85000 )
            {
                th = 0.9;
            }
        }
        if( VOI->GetLargestPossibleRegion().GetSize(2) > 400 )
        {
            if( trachea_voxels > 170000 )
            {
                th = 0.5;
            }
            else if( trachea_voxels > 140000 && trachea_voxels <= 170000 )
            {
                th = 0.75;
            }
            else if( trachea_voxels <= 140000 )
            {
                th = 0.9;
            }
        }
    }
    else if( reconKernel == "LUNG" || reconKernel == "B50f" || reconKernel == "B60f"
             || reconKernel == "D" || reconKernel == "FC50" || reconKernel == "FC52" )
    {
        if( VOI->GetLargestPossibleRegion().GetSize(2) <= 300 )
        {
            if( trachea_voxels > 85000 )
            {
                th = 0.2;
            }
            else if( trachea_voxels > 75000 && trachea_voxels <= 85000 )
            {
                th = 0.3;
            }
            else if( trachea_voxels > 35000 && trachea_voxels <= 75000 )
            {
                th = 0.35;
            }
            else if( trachea_voxels > 10000 && trachea_voxels <= 35000 )//( trachea_voxels > 40000 && trachea_voxels <= 100000 )
            {
                th = 0.5;
            }
            else if( trachea_voxels <= 10000 )
            {
                th = 0.8;
            }
        }
        if( VOI->GetLargestPossibleRegion().GetSize(2) > 300 && VOI->GetLargestPossibleRegion().GetSize(2) <= 400 )
        {
            if( trachea_voxels > 120000 )
            {
                th = 0.2;
            }
            else if( trachea_voxels > 100000 && trachea_voxels <= 120000 ) 
            {
                th = 0.35;
            }
            else if( trachea_voxels > 85000 && trachea_voxels <= 100000 )
            {
                th = 0.5;
            }
            else if( trachea_voxels <= 85000 )
            {
                th = 0.75;
            }
        }
        if( VOI->GetLargestPossibleRegion().GetSize(2) > 400 )
        {
            if( trachea_voxels > 140000 )
            {
                th = 0.2;
            }
            else if( trachea_voxels > 115000 && trachea_voxels <= 140000 )
            {
                th = 0.35;
            }
            else if( trachea_voxels > 80000 && trachea_voxels <= 115000 )
            {
                th = 0.5;
            }
            else if( trachea_voxels <= 80000 )
            {
                th = 0.75;
            }
        }
    }
    else if( reconKernel == "B70f" || reconKernel == "B70s" )
    {
        if(VOI->GetLargestPossibleRegion().GetSize(2) < 300 )
        {
            if( trachea_voxels > 90000)
            {
                th = 0.35;
            }
            else if( trachea_voxels > 60000 && trachea_voxels <= 90000 )
            {
                th = 0.5;
            }
            else if( trachea_voxels > 30000 && trachea_voxels <= 60000 )
            {
                th = 0.6;
            }
            else if( trachea_voxels <= 30000 )
            {
                th = 0.8;
            }
        }
        if(VOI->GetLargestPossibleRegion().GetSize(2) > 300 )
        {
            if( trachea_voxels > 120000)
            {
                th = 0.25;
            }
            else if( trachea_voxels > 80000 && trachea_voxels <= 120000 )
            {
                th = 0.4;
            }
            else if( trachea_voxels > 50000 && trachea_voxels <= 80000 )
            {
                th = 0.6;
            }
            else if( trachea_voxels <= 50000 )
            {
                th = 0.8;
            }
        }
    }

    double g_max         = 1.6; // 0.15 according to Gao's idea
    double g             = 0;

    double n_voxels      = 0;
    double n_voxels_prev = 1;
    double n_voxels_max  = trachea_voxels * th; // To take into account the fact that half trachea is used to mask the input image
        
    /*if( n_voxels_max < 20000 )
    {
        n_voxels_max = 20000;
    }*/

    /** SEGMENTATION PIPELINE */

    typedef itk::ConnectedThresholdImageFilter< InputImageType, InputImageType > ConnectedFilterType; 
    ConnectedFilterType::Pointer thresholdConnected = ConnectedFilterType::New();

    thresholdConnected->SetInput( VOI );			  
    thresholdConnected->SetReplaceValue( labelColor ); 

    InputPixelType UpperThreshold = -930;	

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
    n_voxels = ( StatisticsFilter->GetSum() / labelColor );

    // If the starting threshold gives an empty segmentation, neighbours are checked
    InputImageType::SizeType radius,regionSize;
    InputImageType::IndexType regionIndex;
    InputImageType::RegionType region;	                                                        
    
    if( n_voxels == 0 )
    {
        bool isMinor = 0;
        regionSize.Fill(3);
        regionIndex = index;
        radius.Fill(3);

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
        n_voxels = ( StatisticsFilter->GetSum() / labelColor );
    }

    // If the number of voxels resulting form the segmentation is too high the threshold is iteratively decreased
    if( n_voxels > n_voxels_max )
    {		
        while( n_voxels > n_voxels_max )
        {
            UpperThreshold = UpperThreshold - 10;
            
            thresholdConnected->SetUpper( UpperThreshold );
            thresholdConnected->Update();			
            caster->SetInput( thresholdConnected->GetOutput() );  
            caster->Update();
					
            StatisticsFilter->SetInput(caster->GetOutput());
            StatisticsFilter->Update();			
            n_voxels = ( StatisticsFilter->GetSum() / labelColor );

            if( n_voxels < 5000 )
            {
                bool isMinor = 0;
		
                regionIndex = index;

                regionSize.Fill(3);	
                radius.Fill(3);
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
							
                        StatisticsFilter->SetInput(caster->GetOutput());
                        StatisticsFilter->Update();
                        n_voxels = ( StatisticsFilter->GetSum() / labelColor );
							
                        if( n_voxels > 4000 )
                        {
                            isMinor = 1;
                        }		
                    }
                    counter++;
                }
            }
        }

        // If n_voxels is too small, an increase of the threshold of even 1 HU might cause the violation of g < g_max
        while( g < g_max && n_voxels < n_voxels_max && UpperThreshold <= -800)
        {
            if( n_voxels > 5000 )
            {
                n_voxels_prev = n_voxels;
            }
            UpperThreshold = UpperThreshold + 1;
            thresholdConnected->SetUpper( UpperThreshold );
            thresholdConnected->Update();
            caster->SetInput( thresholdConnected->GetOutput() );  
            caster->Update();
				
            StatisticsFilter->SetInput(caster->GetOutput());
            StatisticsFilter->Update();	
				
            n_voxels = ( StatisticsFilter->GetSum() / labelColor );
            					
            if( n_voxels < 5000 )
            {
                bool isMinor = 0;
		
                regionIndex = index;
                regionSize.Fill(3);	
                radius.Fill(3);
	 
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
				
                        StatisticsFilter->SetInput(caster->GetOutput());
                        StatisticsFilter->Update();
                        n_voxels = ( StatisticsFilter->GetSum() / labelColor );

                        if(n_voxels > 4000)
                        {
                            isMinor = 1;
                        }		
                    }
                    counter++;
                }
            }
            if( n_voxels_prev > 5000 )
            {
                g = double( n_voxels/n_voxels_prev );
            }
        }

        UpperThreshold = UpperThreshold - 1;
        thresholdConnected->SetUpper( UpperThreshold );	
        caster->SetInput( thresholdConnected->GetOutput() );  
        caster->Update();

        StatisticsFilter->SetInput(caster->GetOutput());
        StatisticsFilter->Update();			
        n_voxels = ( StatisticsFilter->GetSum() / labelColor );

        if( n_voxels < 4000 )
        {
            UpperThreshold = UpperThreshold + 1;
            thresholdConnected->SetUpper( UpperThreshold );
            caster->SetInput( thresholdConnected->GetOutput() );  
            caster->Update();
        }
    }
    else      // The threshold is iteratively increased until leakage occurs
    {
        // If n_voxels is too small, an increase of the threshold of even 1 HU might cause the violation of g < g_max
        if( n_voxels_max > 5000 )
        {		
            while( n_voxels < 5000 )
            { 
                UpperThreshold = UpperThreshold + 1;
                thresholdConnected->SetUpper( UpperThreshold );
        
                caster->SetInput( thresholdConnected->GetOutput() );  
                caster->Update();
    			
                StatisticsFilter->SetInput(caster->GetOutput());
                StatisticsFilter->Update();	
    			
                n_voxels = ( StatisticsFilter->GetSum() / labelColor );
            }
        }
        do{	
            UpperThreshold = UpperThreshold + 20;

            n_voxels_prev = n_voxels;	
    			
            thresholdConnected->SetUpper( UpperThreshold );
    	
            caster->SetInput( thresholdConnected->GetOutput() );  
            caster->Update();
    			
            StatisticsFilter->SetInput(caster->GetOutput());
            StatisticsFilter->Update();
    
            n_voxels = ( StatisticsFilter->GetSum() / labelColor );
            g = double( n_voxels/n_voxels_prev );	// double((n_voxels - n_voxels_prev)/n_voxels_prev) according to Gao et al.
        }while( g < g_max && n_voxels < n_voxels_max && UpperThreshold <= -800 );
		
        UpperThreshold = UpperThreshold - 20;	
        thresholdConnected->SetUpper( UpperThreshold );

        caster->SetInput( thresholdConnected->GetOutput() );  
        caster->Update();
	
        StatisticsFilter->SetInput(caster->GetOutput());
        StatisticsFilter->Update();	
        n_voxels = ( StatisticsFilter->GetSum() / labelColor );

        do{	
            UpperThreshold = UpperThreshold + 1;
            n_voxels_prev = n_voxels;	
    			
            thresholdConnected->SetUpper( UpperThreshold );		
    	
            caster->SetInput( thresholdConnected->GetOutput() );  
            caster->Update();
    
            StatisticsFilter->SetInput(caster->GetOutput());
            StatisticsFilter->Update();
			
            n_voxels = ( StatisticsFilter->GetSum() / labelColor );
            g = double( n_voxels/n_voxels_prev );	// double((n_voxels - n_voxels_prev)/n_voxels_prev) according to Gao et al.
        }while( g < g_max && n_voxels < n_voxels_max && UpperThreshold <= -800 );
			
        UpperThreshold = UpperThreshold - 1;

        thresholdConnected->SetUpper( UpperThreshold );	
        caster->SetInput( thresholdConnected->GetOutput() );  
        caster->Update();
    
        StatisticsFilter->SetInput(caster->GetOutput());
        StatisticsFilter->Update();			
        n_voxels = ( StatisticsFilter->GetSum() / labelColor );
    
        if( n_voxels_max > 4000 && n_voxels < 4000 )
        {
            UpperThreshold = UpperThreshold + 1;

            thresholdConnected->SetUpper( UpperThreshold );
            caster->SetInput( thresholdConnected->GetOutput() );  
            caster->Update();
        }
    }
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
        std::cerr << "Exception caught in pasteImage!" << std::endl;
        std::cerr << excep << std::endl;
    } 
		
    return pasteFilter->GetOutput();
}


int main( int argc, char *argv[] )
{

    PARSE_ARGS; 	                         
  	
    typedef itk::ImageFileReader<InputImageType>  ReaderType;				
    ReaderType::Pointer reader = ReaderType::New();     
    
    reader->SetFileName( inputVolume.c_str() );

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

    // Take care of not-supine scanned datasets
    bool flipIm = (reader->GetOutput()->GetDirection()[0][0] == -1 && reader->GetOutput()->GetDirection()[1][1] == -1);

    if( flipIm )
    {
        typedef itk::FlipImageFilter<InputImageType> flipImageFilterType;
        flipImageFilterType::Pointer flipImageFilter = flipImageFilterType::New();
        bool axes[3] = {true,true,false};
        flipImageFilter->SetInput(reader->GetOutput());
        flipImageFilter->SetFlipAxes(axes);
        flipImageFilter->Update();
        reader->GraftOutput(flipImageFilter->GetOutput());
        reader->Update();
    }
	
    // The different labels will be pasted into pasteImage
    OutputImageType::Pointer pasteImage = OutputImageType::New(); 

    pasteImage->SetRegions( reader->GetOutput()->GetRequestedRegion() );                                   
    pasteImage->SetBufferedRegion( reader->GetOutput()->GetBufferedRegion() );
    pasteImage->SetLargestPossibleRegion( reader->GetOutput()->GetLargestPossibleRegion() );
    pasteImage->CopyInformation( reader->GetOutput() );
    pasteImage->Allocate();
    pasteImage->FillBuffer(0);
	
    InputImageType::IndexType tracheaFiducial;
    InputImageType::PointType tracheaPoint;

    std::vector< std::vector<float> > tracheaSeedPoint;
  
    // Finding the trachea seed point 
    int tracheaIndex = 0;

    if( seed.size() == 1 )
    { 		
        tracheaSeedPoint.push_back( seed[tracheaIndex] );

        // Convert to lps the seed point
        tracheaPoint[0] = seed[tracheaIndex][0] * (-reader->GetOutput()->GetDirection()[0][0]);	
        tracheaPoint[1] = seed[tracheaIndex][1] * (-reader->GetOutput()->GetDirection()[1][1]);
        tracheaPoint[2] = seed[tracheaIndex][2] *   reader->GetOutput()->GetDirection()[2][2];

        // Convert the lps physical point to index
        reader->GetOutput()->TransformPhysicalPointToIndex( tracheaPoint, tracheaFiducial );	     		
    }
    else
    {
        if( seed.size() == 0 )
        {
            std::cerr << "No seeds specified!" << std::endl;
            return -1;
        }
        else
        {
            std::cerr << "Please place only one seed point within the trachea!" << std::endl;
            return -1;
        }
    }   
  
    /** TRACHEA SEGMENTATION */
    InputImageType::SizeType cropSize;	                    
    InputImageType::IndexType tracheaCropStart;              

    cropSize[0] = 70; //TO BE FIXED
    tracheaCropStart[0] = tracheaFiducial[0] - 35;   //TO BE FIXED
    if( reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) > 200 )
    {
        cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - 100;
        tracheaCropStart[2] = 100;
    }
    else
    {
        cropSize[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - (reader->GetOutput()->GetLargestPossibleRegion().GetSize(2)/4); //TO BE FIXED
        tracheaCropStart[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2)/4; // TO BE FIXED
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
    
    double fidPs  = double( FiducialSlice[2] );
    double trSz   = double( cropSize[2] );
    double ratio  = fidPs/trSz;

    if(ratio >= 0.85)
    {
        FiducialSlice[2] = cropSize[2]*0.8;
    }

    trachea = TracheaSegmentation( ROIFilter->GetOutput(), FiducialSlice, tracheaSeedPoint, labelValue );

    typedef itk::BinaryBallStructuringElement< OutputImageType::PixelType, DIM > StructuringElementType;
  	
    StructuringElementType structElement;
    StructuringElementType::SizeType radius;
    radius.Fill( 1 );

    structElement.SetRadius( radius );
    structElement.CreateStructuringElement();
	
    typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > CloseType;
    CloseType::Pointer closing = CloseType::New();

    closing->SetInput( trachea );
    closing->SetKernel( structElement );
    closing->SetForegroundValue( labelValue );
    closing->Update();
    trachea = closing->GetOutput();

    /*FINDING THE CARINA IN THE TRACHEA */
    int          value;
    int          xDist;
    int          yDist; 

    unsigned int yFinalMax       = 0;
    unsigned int yFinalMin       = 0;

    unsigned int xCarinaPosition = 0;
    unsigned int yCarinaPosition = 0;
    unsigned int carinaIdx       = 0;

    InputImageType::IndexType carinaIdxPrevSlice;
    carinaIdxPrevSlice.Fill(0);
    InputImageType::IndexType tmpIdx;
    tmpIdx.Fill(0);

    typedef itk::ImageSliceIteratorWithIndex<OutputImageType> SliceIterator;
    SliceIterator sIt(trachea, trachea->GetLargestPossibleRegion());

    sIt.SetFirstDirection(0);
    sIt.SetSecondDirection(1);

    sIt.GoToBegin();
        
    unsigned int limit = FiducialSlice[2] - 5;
 
    while( sIt.GetIndex()[2] < limit )
    {
        unsigned int yMaxPos     = 0;
        unsigned int yMinPos     = trachea->GetLargestPossibleRegion().GetSize(1); 

        unsigned int yMidPos     = 0; 

        unsigned int yCurrentMax = 0;
        unsigned int yCurrentMin = 0;

        unsigned int prevPos     = 0;
        unsigned int prevLine    = 0;
    
        unsigned int prevDiff    = 0;

        while( !sIt.IsAtEndOfSlice() )
        {
            while( !sIt.IsAtEndOfLine() )
            {
                value = sIt.Get();
                if( value == labelValue )
                {
		    double d = sIt.GetIndex()[1] - prevLine;
                    yDist = abs(d);
                    if( prevLine != 0 && yDist <= 5 )
                    {
                        if( sIt.GetIndex()[1] > yMaxPos )
                        {
                            yMaxPos = sIt.GetIndex()[1];
                        }
                        if( sIt.GetIndex()[1] < yMinPos )
                        {
                            yMinPos = sIt.GetIndex()[1];
                        }
                        prevLine = sIt.GetIndex()[1];

                        unsigned int diff = yMaxPos - yMinPos;
                        if( diff > prevDiff)
                        {
                            yCurrentMax = yMaxPos;
                            yCurrentMin = yMinPos;
                        }
                    }
                    else if( prevLine != 0 && yDist >= 5 )
                    {
                        prevDiff = yCurrentMax - yCurrentMin;

                        yMinPos = trachea->GetLargestPossibleRegion().GetSize(1); 
                        yMaxPos = 0; 
                        if( sIt.GetIndex()[1] > yMaxPos )
                        {
                            yMaxPos = sIt.GetIndex()[1];
                        }
                        if( sIt.GetIndex()[1] < yMinPos )
                        {
                            yMinPos = sIt.GetIndex()[1];
                        }
                    }
                    prevLine = sIt.GetIndex()[1];
                }
                ++sIt;
            }
            sIt.NextLine();
        }

        if( yCurrentMax > yCurrentMin )
        {
            yMidPos = yCurrentMin + (yCurrentMax - yCurrentMin) / 2;
        }

        sIt.GoToBeginOfSlice();
        while( !sIt.IsAtEndOfSlice() )
        {
            while( !sIt.IsAtEndOfLine() )
            {
                value = sIt.Get();
                if( value == labelValue )
                {
                    xDist = sIt.GetIndex()[0] - prevPos;
                    if( prevPos != 0 && sIt.GetIndex()[1] == yMidPos 
                        && xDist >= 10 && xDist < 20
                        && sIt.GetIndex()[0] > int(cropSize[0]/3) )
                    {
                        carinaIdx       = sIt.GetIndex()[2];
                        xCarinaPosition = prevPos + (xDist/2);
                        yCarinaPosition = yMidPos;
                        yFinalMax       = yCurrentMax;
                        yFinalMin       = yCurrentMin;
                    }
                    prevPos = sIt.GetIndex()[0];
                }
                ++sIt;
            }
            sIt.NextLine();
        }
           
        carinaIdxPrevSlice[0] = xCarinaPosition; 
        carinaIdxPrevSlice[1] = yCarinaPosition;
        carinaIdxPrevSlice[2] = carinaIdx;

        sIt.NextSlice();
    }

    cropSize[2] -= carinaIdx;
    tracheaCropStart[2] += carinaIdx;
    FiducialSlice[2] = cropSize[2] - ( reader->GetOutput()->GetLargestPossibleRegion().GetSize(2) - tracheaFiducial[2] ); 

    fidPs  = double( FiducialSlice[2] );
    trSz   = double( cropSize[2] );
    ratio  = fidPs/trSz;

    if(ratio >= 0.85)
    {
        FiducialSlice[2] = cropSize[2]*0.8;
    }

    bool exit = 0;
    sIt.GoToBegin();
        
    while( !exit && !sIt.IsAtEnd() )
    {
        while( !exit && !sIt.IsAtEndOfSlice() )
        {
            while( !exit && !sIt.IsAtEndOfLine() )
            {
                value = sIt.Get();
                if( value == labelValue && sIt.GetIndex()[2] >= (carinaIdx + 10 ) && sIt.GetIndex()[2] <= (cropSize[2] - 5) )
                {
                    if( sIt.GetIndex()[0] == 0 || sIt.GetIndex()[0] == cropSize[0] )
                    {
                        cropSize[0] += 10; 
                        tracheaCropStart[0] -= 5;
                        xCarinaPosition += 5;
                        exit = 1;
                    }
                }
                ++sIt;
            }
            sIt.NextLine();
        }
        sIt.NextSlice();
    }

    DesiredRegion.SetSize(  cropSize  );                                                                
    DesiredRegion.SetIndex( tracheaCropStart );                                                          
 
    ROIFilter->SetInput( reader->GetOutput() );						          
    ROIFilter->SetRegionOfInterest( DesiredRegion );					 
    ROIFilter->Update();	

    trachea = TracheaSegmentation( ROIFilter->GetOutput(), FiducialSlice, tracheaSeedPoint, labelValue );

    OutputImageType::IndexType regionIndex = tracheaCropStart;
   
    typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > outputROIFilterType;  
    outputROIFilterType::Pointer extractCarinaSliceFilter = outputROIFilterType::New();	                 
   
    InputImageType::SizeType sliceSize = trachea->GetLargestPossibleRegion().GetSize();
    sliceSize[2] = 1;
    InputImageType::IndexType sliceIndex;
    sliceIndex.Fill(0);
   
    InputImageType::RegionType sliceRegion;                                                  
    sliceRegion.SetSize(  sliceSize  );                                                                
    sliceRegion.SetIndex( sliceIndex ); 

    extractCarinaSliceFilter->SetInput( trachea );	      
    extractCarinaSliceFilter->SetRegionOfInterest( sliceRegion );					 
    extractCarinaSliceFilter->Update();

    SliceIterator singleSliceIt (extractCarinaSliceFilter->GetOutput(), extractCarinaSliceFilter->GetOutput()->GetLargestPossibleRegion());

    singleSliceIt.SetFirstDirection(0);
    singleSliceIt.SetSecondDirection(1);

    singleSliceIt.GoToBegin();

    unsigned int yRightFiducialMaxPos = 0;
    unsigned int yRightFiducialMinPos = trachea->GetLargestPossibleRegion().GetSize(1);
    unsigned int yRightFiducialPos    = 0; 

    unsigned int yLeftFiducialMaxPos  = 0;
    unsigned int yLeftFiducialMinPos  = trachea->GetLargestPossibleRegion().GetSize(1);
    unsigned int yLeftFiducialPos     = 0; 

    while( !singleSliceIt.IsAtEndOfSlice() )
    {
        while( !singleSliceIt.IsAtEndOfLine() )
        {
            value = singleSliceIt.Get();
            if( value == labelValue )
            {
                if( singleSliceIt.GetIndex()[0] > xCarinaPosition )
                {
                    if( singleSliceIt.GetIndex()[1] >= yFinalMin && singleSliceIt.GetIndex()[1] <= yFinalMax  )
                    {
                            if( singleSliceIt.GetIndex()[1] > yLeftFiducialMaxPos )
                            {
                                yLeftFiducialMaxPos = singleSliceIt.GetIndex()[1];                        
                            }
                            if( singleSliceIt.GetIndex()[1] < yLeftFiducialMinPos )
                            {
                                yLeftFiducialMinPos = singleSliceIt.GetIndex()[1];
                            }
                    }
                }
                else if( singleSliceIt.GetIndex()[0] <= xCarinaPosition ) 
                {
                    if( singleSliceIt.GetIndex()[1] >= yFinalMin && singleSliceIt.GetIndex()[1] <= yFinalMax  )
                    {
                        if( singleSliceIt.GetIndex()[1] > yRightFiducialMaxPos )
                        {
                            yRightFiducialMaxPos = singleSliceIt.GetIndex()[1];
                        }
                        if( singleSliceIt.GetIndex()[1] < yRightFiducialMinPos )
                        {
                            yRightFiducialMinPos = singleSliceIt.GetIndex()[1];
                        }                                             
                    }
                }
            }
            ++singleSliceIt;
        }
        singleSliceIt.NextLine();
    }

    if( yRightFiducialMaxPos > yRightFiducialMinPos )
    {
        yRightFiducialPos = yRightFiducialMinPos + (yRightFiducialMaxPos - yRightFiducialMinPos) / 2;
    }

    if( yLeftFiducialMaxPos > yLeftFiducialMinPos )
    {
        yLeftFiducialPos = yLeftFiducialMinPos + (yLeftFiducialMaxPos - yLeftFiducialMinPos) / 2;
    }

    unsigned int xRightFiducialMaxPos = 0;
    unsigned int xRightFiducialMinPos = trachea->GetLargestPossibleRegion().GetSize(1);
    unsigned int xRightFiducialPos    = 0; 

    unsigned int xLeftFiducialMaxPos  = 0;
    unsigned int xLeftFiducialMinPos  = trachea->GetLargestPossibleRegion().GetSize(0);
    unsigned int xLeftFiducialPos     = 0; 

    singleSliceIt.GoToBegin();
    while( !singleSliceIt.IsAtEndOfSlice() )
    {
        while( !singleSliceIt.IsAtEndOfLine() )
        {
            value = singleSliceIt.Get();
            if( value == labelValue )
            {
                if( singleSliceIt.GetIndex()[1] == yRightFiducialPos )
                {
                    if( singleSliceIt.GetIndex()[0] <= xCarinaPosition)
                    {    
                        if( singleSliceIt.GetIndex()[0] > xRightFiducialMaxPos )
                        {
                            xRightFiducialMaxPos = singleSliceIt.GetIndex()[0];
                        }
                        if( singleSliceIt.GetIndex()[0] < xRightFiducialMinPos )
                        {
                            xRightFiducialMinPos = singleSliceIt.GetIndex()[0];
                        }
                    }
                }
                if( singleSliceIt.GetIndex()[1] == yLeftFiducialPos )
                {
                    if( singleSliceIt.GetIndex()[0] > xCarinaPosition )
                    {
                        if( singleSliceIt.GetIndex()[0] > xLeftFiducialMaxPos )
                        {
                            xLeftFiducialMaxPos = singleSliceIt.GetIndex()[0];
                        }
                        if( singleSliceIt.GetIndex()[0] < xLeftFiducialMinPos )
                        {
                            xLeftFiducialMinPos = singleSliceIt.GetIndex()[0];
                        }
                    }
                }   
            }
            ++singleSliceIt;
        }
        singleSliceIt.NextLine();
    }

    if( xRightFiducialMaxPos > xRightFiducialMinPos )
    {
        xRightFiducialPos = xRightFiducialMinPos + (xRightFiducialMaxPos - xRightFiducialMinPos) / 2;
    } 

    if( xLeftFiducialMaxPos > xLeftFiducialMinPos )
    {
        xLeftFiducialPos = xLeftFiducialMinPos + (xLeftFiducialMaxPos - xLeftFiducialMinPos) / 2;
    }

    InputImageType::IndexType rightFiducial, leftFiducial;

    rightFiducial[0] = tracheaCropStart[0] + xRightFiducialPos;
    rightFiducial[1] = yRightFiducialPos;
    rightFiducial[2] = tracheaCropStart[2]; //+ trachea->GetLargestPossibleRegion().GetSize(2) + carinaIdx;

    leftFiducial[0] = tracheaCropStart[0] + xLeftFiducialPos;
    leftFiducial[1] = yLeftFiducialPos;
    leftFiducial[2] = rightFiducial[2];    

    /** RIGHT AND LEFT AIRWAYS SEGMENTATION */
    typedef itk::MaskNegatedImageFilter< InputImageType, OutputImageType, InputImageType > MaskNegatedImageType;
    MaskNegatedImageType::Pointer maskNegFilter = MaskNegatedImageType::New();
    
    typedef itk::StatisticsImageFilter< OutputImageType > StatisticsImageFilterType;		
    StatisticsImageFilterType::Pointer StatisticsFilter = StatisticsImageFilterType::New();	

    StatisticsFilter->SetInput(trachea);
    StatisticsFilter->Update();
    unsigned int numberOfVoxels = ( StatisticsFilter->GetSum() / labelValue );

    // Use half trachea to mask the input image
    typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > ShapeLabelType;
		
    typedef ShapeLabelType::OutputImageType LabelMapType;
    LabelMapType::Pointer labelMap = LabelMapType::New();
		
    typedef itk::ImageDuplicator<OutputImageType> DuplicatorFilterType;

    /** RIGHT AIRWAY SEGMENTATION */
    OutputImageType::Pointer rightHalfTrachea = OutputImageType::New();
            
    DuplicatorFilterType::Pointer rightDuplicatorFilter = DuplicatorFilterType::New();
    rightDuplicatorFilter->SetInputImage(trachea);
    rightDuplicatorFilter->Update();
    rightHalfTrachea = rightDuplicatorFilter->GetOutput();

    SliceIterator rightSIt(rightHalfTrachea, rightHalfTrachea->GetRequestedRegion());

    rightSIt.SetFirstDirection(0);
    rightSIt.SetSecondDirection(1);

    rightSIt.GoToBegin();

    while( !rightSIt.IsAtEnd() )
    {
        unsigned int yMaxPos     = 0;
        unsigned int yMinPos     = rightHalfTrachea->GetLargestPossibleRegion().GetSize(1); 
    
        unsigned int yMidPos     = 0; 

        unsigned int yCurrentMax = 0;
        unsigned int yCurrentMin = 0;

        unsigned int prevLine    = 0;
    
        unsigned int prevDiff    = 0;

        while( !rightSIt.IsAtEndOfSlice() )
        {
            while( !rightSIt.IsAtEndOfLine() )
            {
                value = rightSIt.Get();
                if( value == labelValue )
                {
		    double d = rightSIt.GetIndex()[1] - prevLine;
                    yDist = abs(d);
                    if( prevLine != 0 && yDist <= 5 )
                    {
                        if( rightSIt.GetIndex()[1] > yMaxPos )
                        {
                            yMaxPos = rightSIt.GetIndex()[1];
                        }
                        if( rightSIt.GetIndex()[1] < yMinPos )
                        {
                            yMinPos = rightSIt.GetIndex()[1];
                        }
                        prevLine = rightSIt.GetIndex()[1];

                        unsigned int diff = yMaxPos - yMinPos;
                        if( diff > prevDiff)
                        {
                            yCurrentMax = yMaxPos;
                            yCurrentMin = yMinPos;
                        }
                    }
                    else if( prevLine != 0 && yDist >= 5 )
                    {
                        prevDiff = yCurrentMax - yCurrentMin;

                        yMinPos = rightHalfTrachea->GetLargestPossibleRegion().GetSize(1); 
                        yMaxPos = 0; 
                        if( rightSIt.GetIndex()[1] > yMaxPos )
                        {
                            yMaxPos = rightSIt.GetIndex()[1];
                        }
                        if( rightSIt.GetIndex()[1] < yMinPos )
                        {
                            yMinPos = rightSIt.GetIndex()[1];
                        }
                    }
                    prevLine = rightSIt.GetIndex()[1];
                }
                ++rightSIt;
            }
            rightSIt.NextLine();
        }

        if( yCurrentMax > yCurrentMin )
        {
            yMidPos = yCurrentMin + (yCurrentMax - yCurrentMin) / 2;
        }
    
        unsigned int xMaxPos  = 0;    
        unsigned int xMinPos  = rightHalfTrachea->GetLargestPossibleRegion().GetSize(0);
        unsigned int xMidPos  = 0;
    
        rightSIt.GoToBeginOfSlice();
        while( !rightSIt.IsAtEndOfSlice() )
        {
            while( !rightSIt.IsAtEndOfLine() )
            {
                value = rightSIt.Get();
                if( rightSIt.GetIndex()[1] == yMidPos && value == labelValue )
                {
                    if( rightSIt.GetIndex()[0] > xMaxPos )
                    {
                        xMaxPos = rightSIt.GetIndex()[0];
                    }
                    if( rightSIt.GetIndex()[0] < xMinPos )
                    {
                        xMinPos = rightSIt.GetIndex()[0];
                    }
                }
                ++rightSIt;
            }
            rightSIt.NextLine();
        }
        if( xMaxPos > xMinPos )
        {
            xMidPos = xMinPos + (xMaxPos - xMinPos) / 2;
        }
        rightSIt.GoToBeginOfSlice();
        while( !rightSIt.IsAtEndOfSlice() )
        {
            while( !rightSIt.IsAtEndOfLine() )
            {
                value = rightSIt.Get();
                if( value == labelValue )
                {
                    if( !flipIm && rightSIt.GetIndex()[0] <= xMidPos )
                    {
                        rightSIt.Set(0);
                    }
                    else if( flipIm && rightSIt.GetIndex()[0] >= xMidPos )
                    {
                        rightSIt.Set(0);            
                    }
                }
                ++rightSIt;
            }
            rightSIt.NextLine();
        }
        rightSIt.NextSlice();
    }

    rightSIt.GoToBegin();
    while( !rightSIt.IsAtEndOfSlice() )
    {
        while( !rightSIt.IsAtEndOfLine() )
        {
            value = rightSIt.Get();
            if( value == labelValue )
            {
                rightSIt.Set(0);
            }
            ++rightSIt;
        }
        rightSIt.NextLine();
    }

    OutputImageType::IndexType idx;
    idx = regionIndex;

    pasteImage->FillBuffer(0);
    pasteImage = Paste<OutputImageType>( rightHalfTrachea, idx, pasteImage );

    OutputImageType::Pointer closeImage = OutputImageType::New(); 

    closeImage->SetRegions(reader->GetOutput()->GetRequestedRegion() );                                   
    closeImage->SetBufferedRegion( reader->GetOutput()->GetBufferedRegion() );
    closeImage->SetLargestPossibleRegion( reader->GetOutput()->GetLargestPossibleRegion() );
    closeImage->CopyInformation( reader->GetOutput() );
    closeImage->Allocate();
    closeImage->FillBuffer(labelValue);

    OutputImageType::SizeType sz = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

    if( !flipIm )
    {
        sz[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(0) - (tracheaCropStart[0] + xCarinaPosition ) + 2;
    }
    else
    {
        sz[0] = tracheaCropStart[0] + xCarinaPosition + 2;
    }

    sz[2] = 6;
    OutputImageType::IndexType stIdx; 
    stIdx.Fill(0);   
                       
    OutputImageType::RegionType rg;                                                  
    rg.SetSize(  sz  );
    rg.SetIndex( stIdx );

    typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > outputROIFilterType;
    outputROIFilterType::Pointer cropCloseImageFilter = outputROIFilterType::New();	                 

    cropCloseImageFilter->SetInput( closeImage );	      
    cropCloseImageFilter->SetRegionOfInterest( rg );					 
    cropCloseImageFilter->Update();

    if(!flipIm)
    {
        stIdx[0] = tracheaCropStart[0] + xCarinaPosition - 2;
        stIdx[2] = rightFiducial[2] - 3;
    }
    else
    {
        stIdx[0] = 0;
        stIdx[2] = leftFiducial[2] - 3;
    }	

    pasteImage = Paste<OutputImageType>( cropCloseImageFilter->GetOutput(), stIdx, pasteImage );

    // The half trachea label is used to mask the input image
    maskNegFilter->SetInput( reader->GetOutput() );
    maskNegFilter->SetMaskImage( pasteImage );
    maskNegFilter->Update();	

    OutputImageType::Pointer rightLung = OutputImageType::New();
    std::string reconstructionKernel = reconstructionKernelType.c_str();
    
    if( flipIm )
    {
        rightLung = RightLeftSegmentation( maskNegFilter->GetOutput(), leftFiducial, reconstructionKernel, numberOfVoxels, labelValue ); 
    }
    else
    {
        rightLung = RightLeftSegmentation( maskNegFilter->GetOutput(), rightFiducial, reconstructionKernel, numberOfVoxels, labelValue ); 
    }
				
    ShapeLabelType::Pointer rightLabelConverter = ShapeLabelType::New();
        
    rightLabelConverter->SetInput( rightLung );
    rightLabelConverter->SetInputForegroundValue( labelValue );
    rightLabelConverter->Update();
		
    labelMap = rightLabelConverter->GetOutput();        

    /** LEFT AIRWAY SEGMENTATION */
    OutputImageType::Pointer leftHalfTrachea = OutputImageType::New();

    DuplicatorFilterType::Pointer leftDuplicatorFilter = DuplicatorFilterType::New();
    leftDuplicatorFilter->SetInputImage(trachea);
    leftDuplicatorFilter->Update();
    leftHalfTrachea = leftDuplicatorFilter->GetOutput();

    SliceIterator leftSIt(leftHalfTrachea, leftHalfTrachea->GetRequestedRegion());

    leftSIt.SetFirstDirection(0);
    leftSIt.SetSecondDirection(1);

    leftSIt.GoToBegin();

    while( !leftSIt.IsAtEndOfSlice() )
    {
        while( !leftSIt.IsAtEndOfLine() )
        {
            value = leftSIt.Get();
            if( value == labelValue )
            {
                leftSIt.Set(0);
            }
            ++leftSIt;
        }
        leftSIt.NextLine();
    }

    leftSIt.GoToBegin();
    while( !leftSIt.IsAtEnd() )
    {
        unsigned int yMaxPos     = 0;
        unsigned int yMinPos     = leftHalfTrachea->GetLargestPossibleRegion().GetSize(1); 
    
        unsigned int yMidPos     = 0; 

        unsigned int yCurrentMax = 0;
        unsigned int yCurrentMin = 0;

        unsigned int prevLine    = 0;
    
        unsigned int prevDiff    = 0;

        while( !leftSIt.IsAtEndOfSlice() )
        {
            while( !leftSIt.IsAtEndOfLine() )
            {
                value = leftSIt.Get();
                if( value == labelValue )
                {
		    double d = leftSIt.GetIndex()[1] - prevLine;
                    yDist = abs(d);
                    if( prevLine != 0 && yDist <= 5 )
                    {
                        if( leftSIt.GetIndex()[1] > yMaxPos )
                        {
                            yMaxPos = leftSIt.GetIndex()[1];
                        }
                        if( leftSIt.GetIndex()[1] < yMinPos )
                        {
                            yMinPos = leftSIt.GetIndex()[1];
                        }
                        prevLine = leftSIt.GetIndex()[1];

                        unsigned int diff = yMaxPos - yMinPos;
                        if( diff > prevDiff)
                        {
                            yCurrentMax = yMaxPos;
                            yCurrentMin = yMinPos;
                        }
                    }
                    else if( prevLine != 0 && yDist >= 5 )
                    {
                        prevDiff = yCurrentMax - yCurrentMin;

                        yMinPos = leftHalfTrachea->GetLargestPossibleRegion().GetSize(1); 
                        yMaxPos = 0; 
                        if( leftSIt.GetIndex()[1] > yMaxPos )
                        {
                            yMaxPos = leftSIt.GetIndex()[1];
                        }
                        if( leftSIt.GetIndex()[1] < yMinPos )
                        {
                            yMinPos = leftSIt.GetIndex()[1];
                        }
                    }
                    prevLine = leftSIt.GetIndex()[1];
                }
                ++leftSIt;
            }
            leftSIt.NextLine();
        }

        if( yCurrentMax > yCurrentMin )
        {
            yMidPos = yCurrentMin + (yCurrentMax - yCurrentMin) / 2;
        }
    
        unsigned int xMaxPos  = 0;    
        unsigned int xMinPos  = leftHalfTrachea->GetLargestPossibleRegion().GetSize(0);
        unsigned int xMidPos  = 0;
    
        leftSIt.GoToBeginOfSlice();
        while( !leftSIt.IsAtEndOfSlice() )
        {
            while( !leftSIt.IsAtEndOfLine() )
            {
                value = leftSIt.Get();
                if( leftSIt.GetIndex()[1] == yMidPos && value == labelValue )
                {
                    if( leftSIt.GetIndex()[0] > xMaxPos )
                    {
                        xMaxPos = leftSIt.GetIndex()[0];
                    }
                    if( leftSIt.GetIndex()[0] < xMinPos )
                    {
                        xMinPos = leftSIt.GetIndex()[0];
                    }
                }
                ++leftSIt;
            }
            leftSIt.NextLine();
        }
        if( xMaxPos > xMinPos )
        {
            xMidPos = xMinPos + (xMaxPos - xMinPos) / 2;
        }
        
        leftSIt.GoToBeginOfSlice();
        while( !leftSIt.IsAtEndOfSlice() )
        {
            while( !leftSIt.IsAtEndOfLine() )
            {
                value = leftSIt.Get();
                if( value == labelValue )
                {
                    if( !flipIm && leftSIt.GetIndex()[0] >= xMidPos )
                    {
                        leftSIt.Set(0);
                    }
                    else if( flipIm && leftSIt.GetIndex()[0] <= xMidPos )
                    {
                        leftSIt.Set(0);            
                    }
                }
                ++leftSIt;
            }
            leftSIt.NextLine();
        }
        leftSIt.NextSlice();
    }

    pasteImage->FillBuffer(0);
    pasteImage = Paste<OutputImageType>( leftHalfTrachea, idx, pasteImage );

    if( !flipIm )
    {
        sz[0] = tracheaCropStart[0] + xCarinaPosition + 2;
    }
    else
    {
        sz[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(0) - (tracheaCropStart[0] + xCarinaPosition ) + 2;
    }
    
    stIdx.Fill(0);

    rg.SetSize(  sz  );
    rg.SetIndex( stIdx );

    cropCloseImageFilter->SetInput( closeImage );	      
    cropCloseImageFilter->SetRegionOfInterest( rg );					 
    cropCloseImageFilter->Update();

    if(!flipIm)
    {
        stIdx[0] = 0;
        stIdx[2] = leftFiducial[2] - 3;
    }
    else
    {
        stIdx[0] = tracheaCropStart[0] + xCarinaPosition - 2;
        stIdx[2] = rightFiducial[2] - 3;
    }			
    pasteImage = Paste<OutputImageType>( cropCloseImageFilter->GetOutput(), stIdx, pasteImage );

    maskNegFilter->SetInput( reader->GetOutput() );
    maskNegFilter->SetMaskImage( pasteImage );
    maskNegFilter->Update();

    OutputImageType::Pointer leftLung = OutputImageType::New();
	
    if( flipIm )
    {
        leftLung = RightLeftSegmentation( maskNegFilter->GetOutput(), rightFiducial, reconstructionKernel, numberOfVoxels, labelValue );
    }
    else
    {
        leftLung = RightLeftSegmentation( maskNegFilter->GetOutput(), leftFiducial, reconstructionKernel, numberOfVoxels, labelValue );
    }
				
    ShapeLabelType::Pointer leftLabelConverter = ShapeLabelType::New();

    leftLabelConverter->SetInput( leftLung );
    leftLabelConverter->SetInputForegroundValue( labelValue );
    leftLabelConverter->Update();  

    typedef itk::MergeLabelMapFilter< LabelMapType > MergeFilterType;	  
    MergeFilterType::Pointer mergeFilter = MergeFilterType::New(); 	  
    mergeFilter->SetMethod( MergeFilterType::PACK );
    mergeFilter->SetInput( labelMap );
    mergeFilter->SetInput( 1, leftLabelConverter->GetOutput() );
    mergeFilter->Update();
    labelMap = mergeFilter->GetOutput();

    pasteImage->FillBuffer(0);
    pasteImage = Paste<OutputImageType>( trachea, regionIndex, pasteImage );

    ShapeLabelType::Pointer labelConverter = ShapeLabelType::New();

    labelConverter->SetInput( pasteImage );
    labelConverter->SetInputForegroundValue( labelValue ); 
    labelConverter->Update();
	  
    for( unsigned int i = 0; i < labelConverter->GetOutput()->GetNumberOfLabelObjects(); i++ )
    {
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

    ////////////////////////////////////////////////////////
    /*regionIndex.Fill(0);
    pasteImage->FillBuffer(0);
    pasteImage = Paste<OutputImageType>( leftLung,regionIndex, pasteImage );*/
    ////////////////////////////////////////////////////////
       
    if( flipIm )
    {
        typedef itk::FlipImageFilter<OutputImageType> flipImageFilterType;
        flipImageFilterType::Pointer flipImageFilter = flipImageFilterType::New();
        bool axes[3] = {true,true,false};
        flipImageFilter->SetInput(pasteImage);
        flipImageFilter->SetFlipAxes(axes);
        flipImageFilter->Update();
        pasteImage = flipImageFilter->GetOutput();
    }

    /** CLOSING AND HOLE FILLING TO IMPROVE THE SEGMENTATION RESULT */
    //typedef itk::BinaryBallStructuringElement< OutputImageType::PixelType, DIM > StructuringElementType;
  	
    StructuringElementType newStructElement;
    StructuringElementType::SizeType newRadius;
    newRadius.Fill( 5 );

    newStructElement.SetRadius( newRadius );
    newStructElement.CreateStructuringElement();
	
    //typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > CloseType;
    CloseType::Pointer newClosing = CloseType::New();

    newClosing->SetInput( pasteImage );
    newClosing->SetKernel( newStructElement );
    newClosing->SetForegroundValue( labelValue );
    newClosing->SetSafeBorder( 1 );
    newClosing->Update();
	
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter< OutputImageType > IterativeFillHolesFilterType;
    IterativeFillHolesFilterType::Pointer HoleFilling = IterativeFillHolesFilterType::New();

    OutputImageType::SizeType FillRadius;

    FillRadius.Fill(1);

    HoleFilling->SetInput( newClosing->GetOutput() );
    HoleFilling->SetRadius( FillRadius );
    HoleFilling->SetBackgroundValue( 0 );
    HoleFilling->SetForegroundValue( labelValue );
    HoleFilling->SetMajorityThreshold( 1 );
    HoleFilling->SetMaximumNumberOfIterations( 10 );
    HoleFilling->Update();

    typedef itk::GrayscaleFillholeImageFilter< OutputImageType, OutputImageType > GSFillHolesFilterType;
    GSFillHolesFilterType::Pointer GSHoleFilling = GSFillHolesFilterType::New();

    GSHoleFilling->SetInput( HoleFilling->GetOutput() );
    GSHoleFilling->SetFullyConnected(1);
    GSHoleFilling->Update();

    /** LABEL CREATION */
    typedef  itk::ImageFileWriter<OutputImageType> WriterType;
    WriterType::Pointer labelImage = WriterType::New();

    labelImage->SetFileName( label.c_str() );
    labelImage->SetInput( GSHoleFilling->GetOutput() );
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

