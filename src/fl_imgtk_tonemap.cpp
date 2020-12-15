#ifdef _MSC_VER
#pragma warning(disable : 4018)  
#pragma warning(disable : 4068)  
#pragma warning(disable : 4244)  
#pragma warning(disable : 4996)  
#endif

#include "fl_imgtk.h"
#include "fl_imgtk_minmax.h"

#ifdef USING_OMP
#include <omp.h>
#endif /// of USING_OMP
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// All source code referenced to FreeImage library 3.
// related in tmoColorConvert, tmoReinhard05, tmpDrago03
//
// tmoDrago03, Tone mapping operator (Drago, 2003)
// - Design and implementation by Herv?Drolon (drolon@infonie.fr)
// - It belong to a part of FreeImage 3
//
// tmoReinhard05, Tone mapping operator (Reinhard, 2005)
// - Design and implementation by
//      * Herv?Drolon (drolon@infonie.fr)
//      * Mihail Naydenov (mnaydenov@users.sourceforge.net)
// - It belong to a part of FreeImage 3
//
////////////////////////////////////////////////////////////////////////////////
#define F_EPSILON       1e-06f
#define F_INFINITE      1e+10f

#define F_LOG05         -0.693147f

#define ILLU_CIE_XR     0.640f
#define ILLU_CIE_YR     0.330f
#define ILLU_CIE_XG     0.300f
#define ILLU_CIE_YG     0.600f
#define ILLU_CIE_XB     0.150f
#define ILLU_CIE_YB     0.060f
#define ILLU_CIE_XW     0.3127f
#define ILLU_CIE_YW     0.3290f

#define ILLU_CIE_D      ( ILLU_CIE_XR * \
                          ( ILLU_CIE_YG - ILLU_CIE_YB) + \
                          ILLU_CIE_XG * ( ILLU_CIE_YB - ILLU_CIE_YR ) + \
                          ILLU_CIE_XB * (ILLU_CIE_YR - ILLU_CIE_YG) )

#define ILLU_CIE_C_RD   ( ( 1.0 / ILLU_CIE_YW ) * \
                          ( ILLU_CIE_XW * ( ILLU_CIE_YG - ILLU_CIE_YB ) - \
                            ILLU_CIE_YW * ( ILLU_CIE_XG - ILLU_CIE_XB ) + \
                            ILLU_CIE_XG * ILLU_CIE_YB - ILLU_CIE_XB * ILLU_CIE_YG) )
#define ILLU_CIE_C_GD   ( ( 1.0 / ILLU_CIE_YW ) * \
                          ( ILLU_CIE_XW * ( ILLU_CIE_YB - ILLU_CIE_YR ) - \
                            ILLU_CIE_YW * ( ILLU_CIE_XB - ILLU_CIE_XR ) - \
                            ILLU_CIE_XR * ILLU_CIE_YB + ILLU_CIE_XB * ILLU_CIE_YR ) )
#define ILLU_CIE_C_BD   ( ( 1.0 / ILLU_CIE_YW ) * \
                          ( ILLU_CIE_XW * ( ILLU_CIE_YR - ILLU_CIE_YG ) - \
                            ILLU_CIE_YW * ( ILLU_CIE_XR - ILLU_CIE_XG ) + \
                            ILLU_CIE_XR * ILLU_CIE_YG - ILLU_CIE_XG * ILLU_CIE_YR ) )

////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    unsigned w;
    unsigned h;
    unsigned d;
    float* pixels;
}fl_imgtk_fimg;

////////////////////////////////////////////////////////////////////////////////

/*
** default reference :

static const float RGB2XYZ[3][3] = {
    { 0.41239083F, 0.35758433F, 0.18048081F },
    { 0.21263903F, 0.71516865F, 0.072192319F },
    { 0.019330820F, 0.11919473F, 0.95053220F }
};

static const float XYZ2RGB[3][3] = {
    { 3.2409699F, -1.5373832F, -0.49861079F },
    { -0.96924376F, 1.8759676F, 0.041555084F },
    { 0.055630036F, -0.20397687F, 1.0569715F }
};
*/


// Table of RGB to XYZ (no white balance)
static const double  RGB2XYZ[3][3] =
{
    {
        ILLU_CIE_XR * ILLU_CIE_C_RD / ILLU_CIE_D,
        ILLU_CIE_XG * ILLU_CIE_C_GD / ILLU_CIE_D,
        ILLU_CIE_XB * ILLU_CIE_C_BD / ILLU_CIE_D
    },
    {
        ILLU_CIE_YR * ILLU_CIE_C_RD / ILLU_CIE_D,
        ILLU_CIE_YG * ILLU_CIE_C_GD / ILLU_CIE_D,
        ILLU_CIE_YB * ILLU_CIE_C_BD / ILLU_CIE_D
    },
    {
        (1.0 - ILLU_CIE_XR - ILLU_CIE_YR) * ILLU_CIE_C_RD / ILLU_CIE_D,
        (1.0 - ILLU_CIE_XG - ILLU_CIE_YG) * ILLU_CIE_C_GD / ILLU_CIE_D,
        (1.0 - ILLU_CIE_XB - ILLU_CIE_YB) * ILLU_CIE_C_BD / ILLU_CIE_D
    }
};

// Table of XYZ to RGB (no white balance)
static const double  XYZ2RGB[3][3] =
{
    {
        ( ILLU_CIE_YG - ILLU_CIE_YB - ILLU_CIE_XB * ILLU_CIE_YG + ILLU_CIE_YB * ILLU_CIE_XG )
        / ILLU_CIE_C_RD,
        ( ILLU_CIE_XB - ILLU_CIE_XG - ILLU_CIE_XB * ILLU_CIE_YG + ILLU_CIE_XG * ILLU_CIE_YB )
        / ILLU_CIE_C_RD,
        ( ILLU_CIE_XG * ILLU_CIE_YB - ILLU_CIE_XB * ILLU_CIE_YG )
        / ILLU_CIE_C_RD
    },
    {
        ( ILLU_CIE_YB - ILLU_CIE_YR - ILLU_CIE_YB * ILLU_CIE_XR + ILLU_CIE_YR * ILLU_CIE_XB )
        / ILLU_CIE_C_GD,
        ( ILLU_CIE_XR - ILLU_CIE_XB - ILLU_CIE_XR * ILLU_CIE_YB + ILLU_CIE_XB * ILLU_CIE_YR )
        / ILLU_CIE_C_GD,
        ( ILLU_CIE_XB * ILLU_CIE_YR - ILLU_CIE_XR * ILLU_CIE_YB )
        / ILLU_CIE_C_GD
    },
    {
        ( ILLU_CIE_YR - ILLU_CIE_YG - ILLU_CIE_YR * ILLU_CIE_XG + ILLU_CIE_YG * ILLU_CIE_XR )
        / ILLU_CIE_C_BD,
        ( ILLU_CIE_XG - ILLU_CIE_XR - ILLU_CIE_XG * ILLU_CIE_YR + ILLU_CIE_XR * ILLU_CIE_YG )
        / ILLU_CIE_C_BD,
        ( ILLU_CIE_XR * ILLU_CIE_YG - ILLU_CIE_XG * ILLU_CIE_YR )
        / ILLU_CIE_C_BD
    }
};

////////////////////////////////////////////////////////////////////////////////

static inline double fl_imgtk_biasFunction(const double b, const double x)
{
    // pow(x, log(bias)/log(0.5)
    return pow (x, b);
}

static inline double fl_imgtk_pade_log(const double x)
{
    if( x < 1.0 )
    {
        return ( x * ( 6.0 + x ) / ( 6.0 + 4.0 * x ) );
    }
    else
    if( x < 2.0 )
    {
        return ( x * ( 6.0 + 0.7662 * x ) / ( 5.9897 + 3.7658 * x ) );
    }

    return log( x + 1.0 );
}

////////////////////////////////////////////////////////////////////////////////

fl_imgtk_fimg* fl_imgtk_RGB2F( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) && ( d >= 3 ) )
    {
        fl_imgtk_fimg* newfimg = new fl_imgtk_fimg;
        if ( newfimg != NULL )
        {
            newfimg->w = w;
            newfimg->h = h;
            newfimg->d = d;

            newfimg->pixels = new float[ w * h * d ];
            if ( newfimg->pixels == NULL )
            {
                delete newfimg;
                return NULL;
            }

            #pragma omp parallel for
            for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
            {
                newfimg->pixels[ cnt ] = (float)ptr[ cnt ] / 255.0f;
            }

            return newfimg;
        }
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk_F2RGB( fl_imgtk_fimg* img, bool dranged = false )
{
    if ( img == NULL )
        return NULL;

	float divf = 1.f;
    unsigned w = img->w;
    unsigned h = img->h;
    unsigned d = img->d;

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        uchar* buff = new uchar[ w * h * d ];

        if ( buff != NULL )
        {
			// if dynamic range adjusting, find maximum range.
			if ( dranged == true )
			{
				#pragma omp parallel for reduction(max:divf)
				for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
				{
					divf = MAX( divf, img->pixels[ cnt ] );
				}
				
				#pragma omp parallel for
				for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
				{
					float fconv  = ( img->pixels[ cnt ] / divf ) * 255.0f;
					if ( fconv > 255.f ) fconv = 255.f;
					buff[ cnt ] = (uchar)fconv;
				}
			}
			else
			{
				#pragma omp parallel for
				for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
				{
					float fconv  = img->pixels[ cnt ] * 255.0f;
					if ( fconv > 255.f ) fconv = 255.f;
					buff[ cnt ] = (uchar)fconv;
				}
			}

            return new Fl_RGB_Image( buff, w, h ,d );
        }
    }

    return NULL;
}

bool fl_imgtk_F2RGB_ex( Fl_RGB_Image* dimg, fl_imgtk_fimg* fimg, bool dranged = false )
{
    if ( (dimg == NULL ) || ( fimg == NULL ) )
        return false;

	float divf = 1.f;
    unsigned w = fimg->w;
    unsigned h = fimg->h;
    unsigned d = fimg->d;

    if ( ( dimg->w() != fimg->w ) || ( dimg->h() != fimg->h )
         || ( dimg->d() != fimg->d ) )
    {
        return false;
    }

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        uchar* buff = (uchar*)dimg->data()[0];

        if ( buff != NULL )
        {
			// if dynamic range adjusting, find maximum range.
			if ( dranged == true )
			{
				#pragma omp parallel for reduction(max:divf)
				for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
				{
					divf = MAX( divf, fimg->pixels[ cnt ] );
				}
				
				#pragma omp parallel for
				for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
				{
					float fconv  = ( fimg->pixels[ cnt ] / divf ) * 255.0f;
					if ( fconv > 255.f ) fconv = 255.f;
					buff[ cnt ] = (uchar)fconv;
				}
			}
			else
			{
				#pragma omp parallel for
				for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
				{
					float fconv  = fimg->pixels[ cnt ] * 255.0f;
					if ( fconv > 255.f ) fconv = 255.f;
					buff[ cnt ] = (uchar)fconv;
				}
			}

			return true;
        }
    }

    return false;
}

Fl_RGB_Image* fl_imgtk_clamp_F2RGB( fl_imgtk_fimg* img )
{
    if ( img == NULL )
        return NULL;

    unsigned w = img->w;
    unsigned h = img->h;
    unsigned d = img->d;

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        uchar* buff = new uchar[ w * h * d ];

        if ( buff != NULL )
        {
            unsigned imgsz = w*h*d;
            #pragma omp parallel for
            for( unsigned cnt=0; cnt<imgsz; cnt++ )
            {
                float conv = ( img->pixels[ cnt ] > 1.0f ) ? 1.0f : img->pixels[ cnt ];
                conv *= 255.f;
                conv += 0.5f;
                if ( conv > 255.f )
                    conv = 255.f;
                buff[ cnt ] = (uchar)( conv );
            }

            return new Fl_RGB_Image( buff, w, h ,d );
        }
    }

    return NULL;
}

bool fl_imgtk_clamp_F2RGB_ex( Fl_RGB_Image* dimg, fl_imgtk_fimg* fimg )
{
    if ( ( dimg == NULL ) || ( fimg == NULL ) )
        return false;;

    unsigned w = fimg->w;
    unsigned h = fimg->h;
    unsigned d = fimg->d;

    if ( ( dimg->w() != fimg->w ) || ( dimg->h() != fimg->h )
         || ( dimg->d() != fimg->d ) )
    {
        return false;
    }

    uchar* buff = (uchar*)dimg->data()[0];

    if ( ( w > 0 ) && ( h > 0 ) && ( buff != NULL ) )
    {
        unsigned imgsz = w*h*d;

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            float conv = ( fimg->pixels[ cnt ] > 1.0f ) ? 1.0f : fimg->pixels[ cnt ];
            buff[ cnt ] = (uchar)( conv * 255.0f + 0.5f );
        }

        return true;
    }

    return false;
}


// Convert in-place floating point RGB data to Yxy.
bool fl_imgtk_F2YXY( fl_imgtk_fimg* img )
{
    if ( img != NULL )
    {
        unsigned imgsz = img->w * img->h;
        #pragma omp parellel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            float tf[3] = {0,0,0};
            float* psrc = &img->pixels[ cnt * img->d ];

            for ( unsigned rpt=0; rpt<3; rpt++ )
            {
                tf[rpt] += (float)(RGB2XYZ[rpt][0] * psrc[0]);
                tf[rpt] += (float)(RGB2XYZ[rpt][1] * psrc[1]);
                tf[rpt] += (float)(RGB2XYZ[rpt][2] * psrc[2]);
            }

            float fW = tf[0] + tf[1] + tf[2];
            float fY = tf[1];

            if ( fW > 0.0f )
            {
                psrc[0] = fY;           /// Y
                psrc[1] = tf[0] / fW;   /// x
                psrc[2] = tf[1] / fW;   /// y
            }
            else
            {
                memset( psrc, 0, 3 * sizeof(float) );
            }
        }

        return true;
    }

    return false;
}

// floating RGB to sRGB conversion.
// https://www.w3.org/Graphics/Color/sRGB
fl_imgtk_fimg* fl_imgtk_F2Y( fl_imgtk_fimg* img )
{
    if ( img == NULL )
        return NULL;

    unsigned w = img->w;
    unsigned h = img->h;
    unsigned d = img->d;

    if ( ( w > 0 ) && ( h > 0 ) && ( d >= 3 ) )
    {
        fl_imgtk_fimg* newfimg = new fl_imgtk_fimg;
        if ( newfimg != NULL )
        {
            newfimg->w = w;
            newfimg->h = h;
            newfimg->d = 1;

            newfimg->pixels = new float[ w * h ];
            if ( newfimg->pixels == NULL )
            {
                delete newfimg;
                return NULL;
            }

            #pragma omp parallel for
            for( unsigned cnt=0; cnt<(w*h); cnt++ )
            {
                float* psrc = &img->pixels[ cnt * d ];

                float luma790 = ( ( psrc[0] * 0.2126f ) +
                                  ( psrc[1] * 0.7152f ) +
                                  ( psrc[2] * 0.0722f ) );

                newfimg->pixels[ cnt ] = ( luma790 > 0.0f ) ? luma790 : 0.0f;
            }

            return newfimg;
        }
    }

    return NULL;
}

bool fl_imgtk_YXY2F( fl_imgtk_fimg* img )
{
    if ( img != NULL )
    {
        unsigned imgsz = img->w * img->h;

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            float tf[3] = {0};
            float* psrc = &img->pixels[ cnt * img->d ];

            float fY = psrc[0]; /// Y (red)
            float fX = 0.0f;
            float fZ = 0.0f;
            tf[1] = psrc[1];  /// x (green)
            tf[2] = psrc[2];  /// y (blue)

            if ( ( fY > F_EPSILON ) && ( tf[1] > F_EPSILON ) && ( tf[2] > F_EPSILON ) )
            {
                fX = ( tf[1] * fY ) / tf[2];
                fZ = ( fX / tf[1] ) - fX - fY;
            }
            else
            {
                fX = F_EPSILON;
                fZ = F_EPSILON;
            }

            psrc[0] = fX;
            psrc[1] = fY;
            psrc[2] = fZ;

            tf[0] = 0;
            tf[1] = 0;
            tf[2] = 0;

            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                tf[rpt] += (float)(XYZ2RGB[rpt][0] * psrc[0]);
                tf[rpt] += (float)(XYZ2RGB[rpt][1] * psrc[1]);
                tf[rpt] += (float)(XYZ2RGB[rpt][2] * psrc[2]);
            }

            psrc[0] = tf[0];
            psrc[1] = tf[1];
            psrc[2] = tf[2];

        }

        return true;
    }

    return false;
}

bool fl_imgtk_luminancefromYXY( fl_imgtk_fimg* img, float &maxlumi, float &minlumi, float &worldlumi )
{
    if ( img != NULL )
    {
        unsigned cmax = img->w * img->h;

        double dsum = 0.0;

        #pragma omp parellel for reduction(+:dsum)
        for( unsigned cnt=0; cnt<cmax; cnt++ )
        {
            float fY = MAX( 0.0f, img->pixels[ cnt * img->d ] );

            maxlumi = ( maxlumi < fY ) ? fY : maxlumi;
            minlumi = ( minlumi < fY ) ? minlumi : fY;

            dsum += log( 2.3e-5f + fY );
        }

        double avglog = dsum / ( (double)img->w * (double)img->h );
        worldlumi = (float)exp( avglog );

        return true;
    }

    return false;
}

bool fl_imgtk_luminancefromY( fl_imgtk_fimg* img, float &maxlumi, float &minlumi, float &avglumi, float &avgloglumi )
{
    if ( img == NULL )
        return false;

    if ( img->d != 1 )
        return false;

    unsigned imgsz = img->w * img->h;

    float maxl   = -1e20f;
    float minl   = 1e20f;
    double lsum  = 0.0;
    double llsum = 0.0;

    #pragma omp parallel for reduction(+:lsum,llsum,maxl) reduction(-:minl)
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        float Y = img->pixels[ cnt ];
        maxl    = ( maxl < Y ) ? Y : maxl;
        minl    = ( ( Y > 0.0f ) && ( minl < Y ) ) ? minl : Y;
        lsum   += Y;
        llsum  += log( 2.3e-5f + Y );
    }

    maxlumi    = maxl;
    minlumi    = minl;
    avglumi    = (float)( lsum / (double)imgsz );
    avgloglumi = (float)exp( llsum / (double)imgsz );

    return true;
}

void fl_imgtk_discard_fimg( fl_imgtk_fimg* &img )
{
    if ( img != NULL )
    {
        if ( img->pixels != NULL )
        {
            delete[] img->pixels;
        }

        delete img;
        img = NULL;
    }
}

////////////////////////////////////////////////////////////////////////////////

bool fl_imgtk_tonemappingreinhard( fl_imgtk_fimg* img, fl_imgtk_fimg* Y, float f, float m, float a, float c )
{
    float lumiavg    = 0.0f;
    float lumilogavg = 0.0f;
    float lumimin    = 1.0f;
    float lumimax    = 1.0f;

    float key        = 0.0f;        /// = low key, high-key.


    if ( ( img == NULL ) || ( Y == NULL ) )
        return false;

    if ( Y->d != 1 )
        return false;

    f = exp( -f );

    if ( ( m == 0.0f ) || ( a != 1.0f ) || ( c != 1.0f ) )
    {
        fl_imgtk_luminancefromY( Y, lumimax, lumimin, lumiavg, lumilogavg );

        key = ( log( lumimax ) - lumilogavg )  / ( log( lumimax ) - log( lumimin ) );

        if ( key < 0.0f )
        {
            key = ( log( lumimax ) -  log( lumilogavg ) ) / ( log( lumimax ) - log( lumimin ) );

            if ( key < 0.0f )
            {
                m = 0.3f;
            }
        }
    }

    m = ( m > 0.0f ) ? m : (float)( 0.3 + 0.7 * pow(key, 1.4) );

    float colMax = -1e6f;
    float colMin = +1e6f;

    unsigned imgsz = img->w * img->h;

    if ( ( a == 1.0f ) && ( c == 0.0f ) )
    {
        #pragma omp parallel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            // = pixel luminance
            float pixllumi = Y->pixels[cnt];

            for( unsigned rpt=0; rpt<img->d; rpt++ )
            {
                float* pixel = &img->pixels[ cnt * img->d + rpt ];

                *pixel /= ( *pixel + pow( f * pixllumi, m ) );

                #pragma omp critical
                {
                    colMax = ( *pixel > colMax ) ? *pixel : colMax;
                    colMin = ( *pixel < colMin ) ? *pixel : colMin;
                }
            }
        }
    }
    else
    {
        double chnlavg[3] = {0,0,0};

        if ( ( a != 1.0f ) && ( c != 0.0f ) )
        {
            #pragma omp parallel for
            for( unsigned cnt=0; cnt<imgsz; cnt++ )
            {
                for( unsigned rpt=0; rpt<3; rpt++ )
                {
                    chnlavg[ rpt ] += img->pixels[ cnt * img->d + rpt ];
                }
            }

            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                chnlavg[ rpt ] /= (float)imgsz;
            }
        }

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            float pixllumi = Y->pixels[cnt];

            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                float* pixel  = &img->pixels[ cnt * img->d + rpt ];
                // global light adaptation
                float lightloc = c * *pixel
                                 + ( 1.0f - c ) * pixllumi;

                // local light adaptation
                float lightglb = c * chnlavg[ rpt ] + ( 1.0f - c ) * lumiavg;

                // interpolated pixel light adaptation
                float lightinp = a * lightloc + ( 1.0f - a ) * lightglb;

                *pixel /= *pixel + pow( f * lightinp, m );

                #pragma omp critical
                {
                     colMax = (*pixel > colMax) ? *pixel : colMax;
                     colMin = (*pixel < colMin) ? *pixel : colMin;
                }
            }
        }

    }

    if ( colMax != colMin )
    {
        float colRange = colMax - colMin;

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                float* pixel = &img->pixels[ cnt * img->d + rpt ];

                *pixel = (*pixel - colMin) / colRange;
            }
        }
    }

    return true;
}

bool fl_imgtk_tonemappingdrago( fl_imgtk_fimg* img, float maxLum, float avgLum, float biasParam, float exposure)
{
    if ( img != NULL )
    {
        if ( biasParam == 0.0f )
            biasParam = 0.85f;

        double lmax = maxLum / avgLum;
        double biasp = log( biasParam ) / F_LOG05;
        double divider = log10( lmax + 1.0 );

        unsigned imgsz = img->w * img->h;

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            double Yw = img->pixels[ cnt * img->d ] / avgLum;
            Yw *= exposure;
            double interpol = log( 2.0 + fl_imgtk_biasFunction( biasp, Yw / lmax ) * 8.0 );
            double logv = fl_imgtk_pade_log( Yw );
            img->pixels[ cnt* img->d ] = (float)( (logv / interpol) / divider);
        }

        return true;
    }

    return false;
}

bool fl_imgtk_rec709gammacorrect( fl_imgtk_fimg* img, float gval )
{
    if ( img != NULL )
    {
        float slope = 4.5f;
        float start = 0.018f;
        float fgamma = ( 0.45f / gval ) * 2.0f;

        if ( gval >= 2.1f )
        {
            start = 0.018f / ( gval - 2.0f ) * 7.5f;
            slope = 4.5f * ( gval - 2.0f ) * 7.5f;
        }
        else
        if ( gval <= 1.9f)
        {
            start = 0.018f / ( 2.0f - gval ) * 7.5f;
            slope = 4.5f * ( 2.0f - gval ) * 7.5f;
        }

        unsigned imgsz = img->w * img->h * img->d;

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
			/*
            float* pixel = &img->pixels[cnt];

            *pixel = ( *pixel <= start ) ? *pixel * slope
                     : ( 1.099f * pow( *pixel, fgamma) - 0.099f );
			*/
            float pixel = img->pixels[cnt];

            img->pixels[cnt] = ( pixel <= start ) ? pixel * slope
                               : ( 1.099f * pow( pixel, fgamma) - 0.099f );			
        }

        return true;
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Fl_RGB_Image* fl_imgtk::tonemapping_reinhard( Fl_RGB_Image* src, float intensity, float contrast, float adaptation, float color_correction )
{
    if ( src != NULL )
    {
        fl_imgtk_fimg* fimg = fl_imgtk_RGB2F( src );
        if ( fimg != NULL )
        {
            fl_imgtk_fimg* yimg = fl_imgtk_F2Y( fimg );
            if ( yimg != NULL )
            {
                bool retb = fl_imgtk_tonemappingreinhard( fimg, yimg,
                                                          intensity,
                                                          contrast,
                                                          adaptation,
                                                          color_correction );

                if ( retb == true )
                {
                    Fl_RGB_Image* newimg = fl_imgtk_clamp_F2RGB( fimg );
                    fl_imgtk_discard_fimg( fimg );
                    fl_imgtk_discard_fimg( yimg );

                    return newimg;
                }

                fl_imgtk_discard_fimg( yimg );
            }

            fl_imgtk_discard_fimg( fimg );
        }
    }

    return NULL;
}

bool fl_imgtk::tonemapping_reinhard_ex( Fl_RGB_Image* src, float intensity, float contrast, float adaptation, float color_correction )
{
    if ( src != NULL )
    {
        fl_imgtk_fimg* fimg = fl_imgtk_RGB2F( src );
        if ( fimg != NULL )
        {
            fl_imgtk_fimg* yimg = fl_imgtk_F2Y( fimg );
            if ( yimg != NULL )
            {
                bool retb = fl_imgtk_tonemappingreinhard( fimg, yimg,
                                                          intensity,
                                                          contrast,
                                                          adaptation,
                                                          color_correction );

                if ( retb == true )
                {
                    bool retb = fl_imgtk_clamp_F2RGB_ex( src, fimg );

                    if ( retb == true )
                    {
                        src->uncache();
                    }

                    fl_imgtk_discard_fimg( fimg );
                    fl_imgtk_discard_fimg( yimg );

                    return retb;
                }

                fl_imgtk_discard_fimg( yimg );
            }

            fl_imgtk_discard_fimg( fimg );
        }
    }

    return false;
}


Fl_RGB_Image* fl_imgtk::tonemapping_drago( Fl_RGB_Image* src, float gamma, float exposure )
{
    if ( src == NULL )
        return NULL;

    Fl_RGB_Image* retimg = NULL;

    fl_imgtk_fimg* convimg = fl_imgtk_RGB2F( src );

    if ( convimg != NULL )
    {
        float biasparam = 0.85f;
        float exposparam = (float)pow(2.0, exposure);

        if ( fl_imgtk_F2YXY( convimg ) == true )
        {
            float lumimax = 0.0f;
            float lumimin = 0.0f;
            float lumiavg = 0.0f;

            if ( fl_imgtk_luminancefromYXY( convimg, lumimax, lumimin, lumiavg ) == true )
            {
                if ( fl_imgtk_tonemappingdrago( convimg, lumimax, lumiavg, biasparam, exposparam ) == true )
                {
                    if ( fl_imgtk_YXY2F( convimg ) == true )
                    {
                        if ( gamma != 1.0 )
                        {
                            fl_imgtk_rec709gammacorrect( convimg, gamma );
                        }

                        retimg = fl_imgtk_F2RGB( convimg, false );
                    }
                }
            }
        }

        fl_imgtk_discard_fimg( convimg );
    }

    return retimg;
}

bool fl_imgtk::tonemapping_drago_ex( Fl_RGB_Image* src, float gamma, float exposure )
{
    if ( src == NULL )
        return false;

    fl_imgtk_fimg* convimg = fl_imgtk_RGB2F( src );

    bool retb = false;

    if ( convimg != NULL )
    {
        float biasparam = 0.85f;
        float exposparam = (float)pow(2.0, exposure);

        if ( fl_imgtk_F2YXY( convimg ) == true )
        {
            float lumimax = 0.0f;
            float lumimin = 0.0f;
            float lumiavg = 0.0f;

            if ( fl_imgtk_luminancefromYXY( convimg, lumimax, lumimin, lumiavg ) == true )
            {
                if ( fl_imgtk_tonemappingdrago( convimg, lumimax, lumiavg, biasparam, exposparam ) == true )
                {
                    if ( fl_imgtk_YXY2F( convimg ) == true )
                    {
                        if ( gamma != 1.0 )
                        {
                            fl_imgtk_rec709gammacorrect( convimg, gamma );
                        }

                        retb = fl_imgtk_F2RGB_ex( src, convimg, false );

                        if ( retb == true )
                        {
                            src->uncache();
                        }
                    }
                }
            }
        }

        fl_imgtk_discard_fimg( convimg );
    }

    return retb;
}
