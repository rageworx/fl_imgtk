#include "fl_imgtk.h"
#include "fl_imgtk_minmax.h"

#ifdef USING_OMP
#include <omp.h>
#endif /// of USING_OMP
#include <cmath>

////////////////////////////////////////////////////////////////////////////////

#define CLAHE_MAX_REG_W             256
#define CLAHE_MAX_REG_H             256
#define CLAHE_MAX_RANGE             512

////////////////////////////////////////////////////////////////////////////////
// Some CLAHE depends functions here:

/* CLAHE_ClipHistogram() :
 * This function performs clipping of the histogram and redistribution of bins.
 * The histogram is clipped and the number of excess pixels is counted. Afterwards
 * the excess pixels are equally redistributed across the whole histogram (providing
 * the bin count is smaller than the cliplimit).
 */
void CLAHE_ClipHistogram( unsigned* pHisto, unsigned greyLvl, unsigned clipLimit )
{
    unsigned* pRangePtr = pHisto;
	unsigned* pEndPtr = NULL;
	unsigned* pHistPtr = NULL;

    unsigned excessSz = 0;
	unsigned upperSz;
	unsigned rangeInc;
	unsigned stepsz;
	unsigned cnt;

    long lBinExcess;

	/* calculate total number of excess pixels */
    for ( cnt = 0; cnt < greyLvl; cnt++)
	{
		lBinExcess = (long) pRangePtr[cnt] - (long) clipLimit;

		/* excess in current bin */
		if ( lBinExcess > 0 )
		{
			excessSz += lBinExcess;
		}
    }

    /* Second part: clip histogram and redistribute excess pixels in each bin */
    rangeInc = excessSz / greyLvl;		 /// average binincrement
    upperSz  = clipLimit - rangeInc;	 /// Bins larger than upperSz set to cliplimit

    for ( cnt=0; cnt<greyLvl; cnt++ )
	{
		if (pHisto[cnt] > clipLimit)
		{
			pHisto[cnt] = clipLimit;    /// clip bin
		}
		else
		{
			if (pHisto[cnt] > upperSz)
			{	/* high bin count */
				excessSz    -= pHisto[cnt] - upperSz;
				pHisto[cnt]  = clipLimit;
			}
			else
			{	/* low bin count */
				excessSz    -= rangeInc;
				pHisto[cnt] += rangeInc;
			}
		}
    }

    while ( excessSz > 0 )
	{   /* Redistribute remaining excess  */
		pEndPtr = &pHisto[greyLvl]; pHistPtr = pHisto;

		while ( ( excessSz > 0 ) && ( pHistPtr < pEndPtr )  )
		{
			stepsz = greyLvl / excessSz;

			if ( stepsz < 1 )
			{
				/* stepsize at least 1 */
				stepsz = 1;
			}

			for ( pRangePtr=pHistPtr;
			      pRangePtr < pEndPtr && excessSz;
				  pRangePtr += stepsz)
			{
				if (*pRangePtr < clipLimit)
				{
					/* reduce excess */
					(*pRangePtr)++;
					excessSz--;
				}
			}

			/* restart redistributing on other bin location */
			pHistPtr++;
		}
    }
}

/* CLAHE_MakeHistogram() :
 * This function classifies the greylevels present in the array image into
 * a greylevel histogram. The pLookupTable specifies the relationship
 * between the greyvalue of the pixel (typically between 0 and 4095) and
 * the corresponding bin in the histogram (usually containing only 128 bins).
 */
void CLAHE_MakeHistogram ( uchar* pImage,
                           unsigned imgWidth,
                           unsigned rgnszW, unsigned rgnszH,
                           unsigned* pHisto,
                           unsigned greyLvl, uchar* pLUT )
{
    uchar*   pImgPtr = NULL;
    unsigned cnt;

	/* clear histogram */
    for ( cnt=0; cnt<greyLvl; cnt++ )
	{
		pHisto[cnt] = 0L;
	}

    for ( cnt=0; cnt< rgnszH; cnt++ )
	{
		pImgPtr = &pImage[rgnszW];

		while ( pImage < pImgPtr )
		{
			pHisto[ pLUT[ *pImage++ ] ]++;
		}

		pImgPtr += imgWidth;
		pImage = &pImgPtr[-(long)rgnszW];
    }
}

/* CLAHE_MapHistogram() :
 * This function calculates the equalized lookup table (mapping) by
 * cumulating the input histogram. Note: lookup table is rescaled in range [Min..Max].
 */
void CLAHE_MapHistogram ( unsigned* pHisto,
                          uchar Min, uchar Max,
                          unsigned greyLvl, unsigned pixelsz )
{
    unsigned cnt = 0;
	unsigned sum = 0;

    const float    fScale = ( (float)(Max - Min) ) / pixelsz;
    const unsigned ulMin  = (unsigned) Min;

    for ( cnt=0; cnt<greyLvl; cnt++)
	{
		sum += pHisto[cnt];
		pHisto[cnt] = (unsigned)( ulMin + sum * fScale );

		if ( pHisto[cnt] > Max )
		{
			pHisto[cnt] = Max;
		}
    }
}

/* CLAHE_MakeLut() :
 * To speed up histogram clipping, the input image [Min,Max] is scaled down to
 * [0,uiNrBins-1]. This function calculates the LUT.
 */
void CLAHE_MakeLut( uchar* pLUT, uchar Min, uchar Max, unsigned ranges )
{
    const uchar BinSize = (uchar) (1 + (Max - Min) / ranges);

    for ( unsigned cnt=Min; cnt<=Max; cnt++)
    {
        pLUT[cnt] = (cnt - Min) / BinSize;
    }
}

/* CLAHE_Interpolate() :
 * pImage      - pointer to input/output image
 * imgWidth      - resolution of image in x-direction
 * pMap*     - mappings of greylevels from histograms
 * subszW     - subszW of image submatrix
 * subszH     - subszH of image submatrix
 * pLUT	       - lookup table containing mapping greyvalues to bins
 * This function calculates the new greylevel assignments of pixels within a submatrix
 * of the image with size subszW and subszH. This is done by a bilinear interpolation
 * between four different mappings in order to eliminate boundary artifacts.
 * It uses a division; since division is often an expensive operation, I added code to
 * perform a logical shift instead when feasible.
 */
void CLAHE_Interpolate( uchar* pImage,
                        int imgWidth, unsigned* pMapLU,
                        unsigned* pMapRU, unsigned* pMapLB, unsigned* pMapRB,
                        unsigned subszW, unsigned subszH, uchar* pLUT)
{
    const unsigned incSz = imgWidth-subszW; /* Pointer increment after processing row */
    uchar GreyValue;
	unsigned normFactor = subszW * subszH; /* Normalization factor */

    unsigned coefW = 0;
	unsigned coefH = 0;
	unsigned invcoefW = 0;
	unsigned invcoefH = 0;
	unsigned shifts = 0;

	/* If normFactor is not a power of two, use division */
    if ( normFactor & (normFactor - 1) )
	{
		for ( coefH=0, invcoefH = subszH;
		      coefH < subszH;
		      coefH++, invcoefH--,pImage+=incSz )
		{
			for ( coefW=0, invcoefW = subszW;
			      coefW < subszW;
			      coefW++, invcoefW--)
			{
				/* get histogram bin value */
				GreyValue = pLUT[*pImage];

				*pImage++ = (uchar)( ( invcoefH *
									 ( invcoefW*pMapLU[GreyValue] +  coefW * pMapRU[GreyValue] ) +
						               coefH * (invcoefW * pMapLB[GreyValue] +
                                       coefW * pMapRB[GreyValue])) / normFactor );
			}
		}
	}
    else
	{	/* avoid the division and use a right shift instead */
		while ( normFactor >>= 1 )
		{
			/* Calculate 2log of normFactor */
			shifts++;
		}

		for ( coefH = 0, invcoefH = subszH;
		      coefH < subszH;
			  coefH++, invcoefH--,pImage+=incSz )
		{
			 for ( coefW = 0, invcoefW = subszW;
			       coefW < subszW;
			       coefW++, invcoefW-- )
			{
				/* get histogram bin value */
				GreyValue = pLUT[*pImage];

				*pImage++ = (uchar)( ( invcoefH *
				                       ( invcoefW * pMapLU[GreyValue] + coefW * pMapRU[GreyValue] ) +
                                         coefH * (invcoefW * pMapLB[GreyValue] +
										 coefW * pMapRB[GreyValue])) >> shifts );
			}
		}
    }
}

////////////////////////////////////////////////////////////////////////////////

/* applyCLAHE() :
 *   pImage - Pointer to the input/output image
 *   imgWidth - Image resolution in the X direction
 *   imgHeight - Image resolution in the Y direction
 *   Min - Minimum greyvalue of input image (also becomes minimum of output image)
 *   Max - Maximum greyvalue of input image (also becomes maximum of output image)
 *   rgnWidth - Number of contextial regions in the X direction (min 2, max DEF_MAX_WIDTH)
 *   rgnHeight - Number of contextial regions in the Y direction (min 2, max DEF_MAX_HEIGHT)
 *   ranges - Number of greybins for histogram ("dynamic range")
 *   float fCliplimit - Normalized cliplimit (higher values give more contrast)
 * The number of "effective" greylevels in the output image is set by ranges; selecting
 * a small value (eg. 128) speeds up processing and still produce an output image of
 * good quality. The output image will have the same minimum and maximum value as the input
 * image. A clip limit smaller than 1 results in standard (non-contrast limited) AHE.
 */

bool applyCLAHE( uchar* pImage,
                 unsigned imgWidth, unsigned imgHeight,
                 uchar Min, uchar Max,
                 unsigned rgnWidth, unsigned rgnHeight,
                 unsigned ranges, float fCliplimit )
{
    /* counters */
    unsigned cnt_x;
    unsigned cnt_y;

    /* size of context. reg. and subimages */
    unsigned subszW = 0;
    unsigned subszH = 0;
    unsigned subImgW = 0;
    unsigned subImgH = 0;

    /* auxiliary variables interpolation routine */
    unsigned cnt_xL = 0;
    unsigned cnt_xR = 0;
    unsigned cnt_yU = 0;
    unsigned cnt_yB = 0;

    /* clip limit and region pixel count */
    unsigned clipLimit = 0;
    unsigned pixelCnts = 0;

    /* pointer to image */
    uchar* pImgPtr = NULL;

    /* lookup table used for scaling of input image */
    uchar aLUT[ 256 ] = {0};

    /* pointer to histogram and mappings*/
    unsigned* pulHist = NULL;
    unsigned* pMapArray = NULL;

    /* auxiliary pointers interpolation */
    unsigned* pulLU = NULL;
    unsigned* pulLB = NULL;
    unsigned* pulRU = NULL;
    unsigned* pulRB = NULL;

    if ( rgnWidth > CLAHE_MAX_REG_W )
		rgnWidth = CLAHE_MAX_REG_W;

    if ( rgnHeight > CLAHE_MAX_REG_H )
		rgnHeight = CLAHE_MAX_REG_H;

    if ( ( imgWidth % rgnWidth ) > 0 )
		return false;       /// x-resolution no multiple of rgnWidth

    if ( ( imgHeight % rgnHeight ) > 0 )
		return false;       /// y-resolution no multiple of rgnHeight

    if ( Min >= Max )
		return false;       /// minimum equal or larger than maximum

    if ( rgnWidth < 2 )
        rgnWidth = 2;

    if ( rgnHeight < 2 )
		rgnHeight = 2;

    if ( fCliplimit == 1.0f )
		return true;	    /// is OK, immediately returns original image.

    if ( ( fCliplimit > 0.0f ) && ( fCliplimit < 1.1f ) )
        fCliplimit = 1.1f;

    if ( ranges == 0 )
		ranges = CLAHE_MAX_RANGE;   /// default value when not specified

    pMapArray = new unsigned[ rgnWidth * rgnHeight * ranges ];

    if ( pMapArray == NULL )
		return false;

    subszW    = imgWidth/rgnWidth; subszH = imgHeight/rgnHeight;  /* Actual size of contextual regions */
    pixelCnts = (unsigned)subszW * (unsigned)subszH;

    if(fCliplimit > 0.0)
	{
		/* Calculate actual cliplimit	 */
		clipLimit = (unsigned) (fCliplimit * (subszW * subszH) / ranges);
        clipLimit = (clipLimit < 1UL) ? 1UL : clipLimit;
    }
    else
	{
		/* Large value, do not clip (AHE) */
		clipLimit = 1UL<<14;
	}

    CLAHE_MakeLut(aLUT, Min, Max, ranges);	  /* Make lookup table for mapping of greyvalues */

    /* Calculate greylevel mappings for each contextual region */
    for ( cnt_y=0, pImgPtr=pImage; cnt_y<rgnHeight; cnt_y++ )
	{
		for ( cnt_x=0; cnt_x<rgnWidth; cnt_x++, pImgPtr+=subszW )
		{
			pulHist = &pMapArray[ranges * (cnt_y * rgnWidth + cnt_x)];

			CLAHE_MakeHistogram( pImgPtr,
                                 imgWidth, subszW, subszH,
                                 pulHist, ranges, aLUT );
			CLAHE_ClipHistogram( pulHist, ranges, clipLimit);
			CLAHE_MapHistogram( pulHist, Min, Max, ranges, pixelCnts);
		}

		pImgPtr += (subszH - 1) * imgWidth;		  /* skip lines, set pointer */
    }

    /* Interpolate greylevel mappings to get CLAHE image */
    for (pImgPtr = pImage, cnt_y = 0; cnt_y <= rgnHeight; cnt_y++)
	{
		if (cnt_y == 0)
		{					  /* special case: top row */
			subImgH = subszH >> 1;  cnt_yU = 0; cnt_yB = 0;
		}
		else
		{
			if (cnt_y == rgnHeight)
			{				  /* special case: bottom row */
				subImgH = subszH >> 1;	cnt_yU = rgnHeight-1;	 cnt_yB = cnt_yU;
			}
			else
			{					  /* default values */
				subImgH = subszH; cnt_yU = cnt_y - 1; cnt_yB = cnt_yU + 1;
			}
		}
		for (cnt_x = 0; cnt_x <= rgnWidth; cnt_x++)
		{
			if (cnt_x == 0)
			{				  /* special case: left column */
				subImgW = subszW >> 1; cnt_xL = 0; cnt_xR = 0;
			}
			else
			{
				if (cnt_x == rgnWidth)
				{			  /* special case: right column */
					subImgW = subszW >> 1;  cnt_xL = rgnWidth - 1; cnt_xR = cnt_xL;
				}
				else
				{					  /* default values */
					subImgW = subszW; cnt_xL = cnt_x - 1; cnt_xR = cnt_xL + 1;
				}
			}

			pulLU = &pMapArray[ranges * (cnt_yU * rgnWidth + cnt_xL)];
			pulRU = &pMapArray[ranges * (cnt_yU * rgnWidth + cnt_xR)];
			pulLB = &pMapArray[ranges * (cnt_yB * rgnWidth + cnt_xL)];
			pulRB = &pMapArray[ranges * (cnt_yB * rgnWidth + cnt_xR)];

			CLAHE_Interpolate( pImgPtr,
                               imgWidth, pulLU, pulRU, pulLB, pulRB,
                               subImgW, subImgH, aLUT );

			pImgPtr += subImgW;			  /* set pointer on next matrix */
		}

		pImgPtr += (subImgH - 1) * imgWidth;
    }

    delete[] pMapArray;

	return true;
}

////////////////////////////////////////////////////////////////////////////////

Fl_RGB_Image* fl_imgtk::CLAHE( Fl_RGB_Image* src, unsigned regionW, unsigned regionH, float cliplimit )
{
    if ( src == NULL )
        return NULL;

    if ( src->d() < 3 )
        return NULL;

    unsigned imgWidth = src->w();
    unsigned imgHeight = src->h();

    if ( ( imgWidth == 0 ) || ( imgHeight == 0 ) )
        return NULL;

    unsigned imgsz = imgWidth * imgHeight;

    uchar* rbuff = (uchar*)src->data()[0];

    // Split R/G/B channels .
    uchar* data_r = new uchar[ imgsz ];
    uchar* data_g = new uchar[ imgsz ];
    uchar* data_b = new uchar[ imgsz ];

    if ( data_r == NULL )
        return NULL;

    if ( data_g == NULL )
    {
        delete[] data_r;
        return NULL;
    }

    if ( data_b == NULL )
    {
        delete[] data_r;
        delete[] data_g;
        return NULL;
    }

    uchar min_rgb[3] = {255,255,255};
    uchar max_rgb[3] = {0,0,0};

    #pragma omp parellel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        uchar* ptr = &rbuff[ cnt * src->d() ];
        data_r[ cnt ] = ptr[ 0 ];
        data_g[ cnt ] = ptr[ 1 ];
        data_b[ cnt ] = ptr[ 2 ];

        if ( min_rgb[ 0 ] > data_r[ cnt ] ) min_rgb[ 0 ] = data_r[ cnt ];
        else
        if ( max_rgb[ 0 ] < data_r[ cnt ] ) max_rgb[ 0 ] = data_r[ cnt ];

        if ( min_rgb[ 1 ] > data_g[ cnt ] ) min_rgb[ 1 ] = data_g[ cnt ];
        else
        if ( max_rgb[ 1 ] < data_g[ cnt ] ) max_rgb[ 1 ] = data_g[ cnt ];

        if ( min_rgb[ 2 ] > data_b[ cnt ] ) min_rgb[ 2 ] = data_b[ cnt ];
        else
        if ( max_rgb[ 2 ] < data_b[ cnt ] ) max_rgb[ 2 ] = data_b[ cnt ];

    }

    // RED->GREEN->BLUE
    // Skip failure.
    applyCLAHE( data_r, imgWidth, imgHeight,
                min_rgb[0], max_rgb[0], regionW, regionH, 256, cliplimit );
    applyCLAHE( data_g, imgWidth, imgHeight,
                min_rgb[1], max_rgb[1], regionW, regionH, 256, cliplimit );
    applyCLAHE( data_b, imgWidth, imgHeight,
                min_rgb[2], max_rgb[2], regionW, regionH, 256, cliplimit );

    Fl_RGB_Image* newimg = (Fl_RGB_Image*)src->copy();

    if ( newimg != NULL )
    {
        rbuff = (uchar*)newimg->data()[0];

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            uchar* ptr = &rbuff[ cnt * newimg->d() ];
            ptr[ 0 ] = data_r[ cnt ];
            ptr[ 1 ] = data_g[ cnt ];
            ptr[ 2 ] = data_b[ cnt ];
        }

        newimg->uncache();
    }

    delete[] data_r;
    delete[] data_g;
    delete[] data_b;

    return newimg;
}

bool fl_imgtk::CLAHE_ex( Fl_RGB_Image* src, unsigned regionW, unsigned regionH, float cliplimit )
{
    if ( src == NULL )
        return false;

    if ( src->d() < 3 )
        return false;

    unsigned imgWidth = src->w();
    unsigned imgHeight = src->h();

    if ( ( imgWidth == 0 ) || ( imgHeight == 0 ) )
        return false;

    unsigned imgsz = imgWidth * imgHeight;

    uchar* rbuff = (uchar*)src->data()[0];

    // Split R/G/B channels .
    uchar* data_r = new uchar[ imgsz ];
    uchar* data_g = new uchar[ imgsz ];
    uchar* data_b = new uchar[ imgsz ];

    if ( data_r == NULL )
        return false;

    if ( data_g == NULL )
    {
        delete[] data_r;
        return false;
    }

    if ( data_b == NULL )
    {
        delete[] data_r;
        delete[] data_g;
        return false;
    }

    uchar min_rgb[3] = {255,255,255};
    uchar max_rgb[3] = {0,0,0};

    #pragma omp parellel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        uchar* ptr = &rbuff[ cnt * src->d() ];
        data_r[ cnt ] = ptr[ 0 ];
        data_g[ cnt ] = ptr[ 1 ];
        data_b[ cnt ] = ptr[ 2 ];

        if ( min_rgb[ 0 ] > data_r[ cnt ] ) min_rgb[ 0 ] = data_r[ cnt ];
        else
        if ( max_rgb[ 0 ] < data_r[ cnt ] ) max_rgb[ 0 ] = data_r[ cnt ];

        if ( min_rgb[ 1 ] > data_g[ cnt ] ) min_rgb[ 1 ] = data_g[ cnt ];
        else
        if ( max_rgb[ 1 ] < data_g[ cnt ] ) max_rgb[ 1 ] = data_g[ cnt ];

        if ( min_rgb[ 2 ] > data_b[ cnt ] ) min_rgb[ 2 ] = data_b[ cnt ];
        else
        if ( max_rgb[ 2 ] < data_b[ cnt ] ) max_rgb[ 2 ] = data_b[ cnt ];

    }

    // RED->GREEN->BLUE
    // Skip failure.
    applyCLAHE( data_r, imgWidth, imgHeight,
                min_rgb[0], max_rgb[0], regionW, regionH, 256, cliplimit );
    applyCLAHE( data_g, imgWidth, imgHeight,
                min_rgb[1], max_rgb[1], regionW, regionH, 256, cliplimit );
    applyCLAHE( data_b, imgWidth, imgHeight,
                min_rgb[2], max_rgb[2], regionW, regionH, 256, cliplimit );

	#pragma omp parellel for
	for( unsigned cnt=0; cnt<imgsz; cnt++ )
	{
		uchar* ptr = &rbuff[ cnt * src->d() ];
		ptr[ 0 ] = data_r[ cnt ];
		ptr[ 1 ] = data_g[ cnt ];
		ptr[ 2 ] = data_b[ cnt ];
	}

    delete[] data_r;
    delete[] data_g;
    delete[] data_b;

	src->uncache();

    return true;
}

Fl_RGB_Image* fl_imgtk::noire( Fl_RGB_Image* src, unsigned regionW, unsigned regionH, float cliplimit, float bright )
{
    if ( src == NULL )
        return NULL;

    if ( src->d() < 3 )
        return NULL;

    unsigned imgWidth = src->w();
    unsigned imgHeight = src->h();

    if ( ( imgWidth == 0 ) || ( imgHeight == 0 ) )
        return NULL;

    if ( bright < 0.1f )
        bright = 0.1f;

    if ( cliplimit < 1.0f )
        cliplimit = 1.1f;

    unsigned imgsz = imgWidth * imgHeight;

    uchar* rbuff = (uchar*)src->data()[0];

    // RGB average
    uchar* data_rgb_avr = new uchar[ imgsz ];

    if ( data_rgb_avr == NULL )
        return NULL;

    uchar min_rgb = 255;
    uchar max_rgb = 0;

    #pragma omp parellel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        uchar* ptr = &rbuff[ cnt * src->d() ];

        float rgb_avr = (float)( ptr[ 0 ] + ptr[ 1 ] + ptr[ 2 ] ) / 3.0f;

        data_rgb_avr[ cnt ] = (unsigned)rgb_avr;

        if ( (float)min_rgb > rgb_avr ) min_rgb = (uchar)rgb_avr;
        else
        if ( (float)max_rgb < rgb_avr ) max_rgb = (uchar)rgb_avr;
    }

    bool retb = applyCLAHE( data_rgb_avr, imgWidth, imgHeight,
                            min_rgb, max_rgb, regionW, regionH,
                            255, cliplimit );

    if ( retb == false )
    {
        delete[] data_rgb_avr;
        return NULL;
    }

    Fl_RGB_Image* newimg = (Fl_RGB_Image*)src->copy();

    if ( newimg != NULL )
    {
        rbuff = (uchar*)newimg->data()[0];

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            uchar* ptr = &rbuff[ cnt * newimg->d() ];

            float lumif = ( (float)data_rgb_avr[cnt] / 255.0f ) * bright;

            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                ptr[ rpt ] = (uchar) MIN( 255.0f,  (float)ptr[ rpt ] * lumif + 0.5f );
            }
        }
    }

    delete[] data_rgb_avr;

    return newimg;
}

bool fl_imgtk::noire_ex( Fl_RGB_Image* src, unsigned regionW, unsigned regionH, float cliplimit, float bright )
{
    if ( src == NULL )
        return false;

    if ( src->d() < 3 )
        return false;

    unsigned imgWidth = src->w();
    unsigned imgHeight = src->h();

    if ( ( imgWidth == 0 ) || ( imgHeight == 0 ) )
        return false;

    if ( bright < 0.1f )
        bright = 0.1f;

    if ( cliplimit < 1.0f )
        cliplimit = 1.1f;

    unsigned imgsz = imgWidth * imgHeight;

    uchar* rbuff = (uchar*)src->data()[0];

    // RGB average
    uchar* data_rgb_avr = new uchar[ imgsz ];

    if ( data_rgb_avr == NULL )
        return NULL;

    uchar min_rgb = 255;
    uchar max_rgb = 0;

    #pragma omp parellel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        uchar* ptr = &rbuff[ cnt * src->d() ];

        float rgb_avr = (float)( ptr[ 0 ] + ptr[ 1 ] + ptr[ 2 ] ) / 3.0f;

        data_rgb_avr[ cnt ] = (unsigned)rgb_avr;

        if ( (float)min_rgb > rgb_avr ) min_rgb = (uchar)rgb_avr;
        else
        if ( (float)max_rgb < rgb_avr ) max_rgb = (uchar)rgb_avr;
    }

    bool retb = applyCLAHE( data_rgb_avr, imgWidth, imgHeight,
                            min_rgb, max_rgb, regionW, regionH,
                            255, cliplimit );

    if ( retb == false )
    {
        delete[] data_rgb_avr;
        return false;
    }

    #pragma omp parellel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        uchar* ptr = &rbuff[ cnt * src->d() ];

        float lumif = ( (float)data_rgb_avr[cnt] / 255.0f ) * bright;

        for( unsigned rpt=0; rpt<3; rpt++ )
        {
            ptr[ rpt ] = (uchar) MIN( 255.0f,  (float)ptr[ rpt ] * lumif + 0.5f );
        }
    }

    delete[] data_rgb_avr;

    return true;
}

