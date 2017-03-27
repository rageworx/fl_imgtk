#include "fl_imgtk.h"
#include "fl_imgtk_minmax.h"

#include <omp.h>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// All source code referenced to FreeImage library 3.
// related in tmoColorConvert, tmoReinhard05, tmpDrago03

#define F_EPSILON       1e-06f
#define F_INFINITE      1e+10f

#define F_LOG05         -0.693147

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

#define ILLU_CIE_C_RD   ( ( 1.0f / ILLU_CIE_YW ) * \
                          ( ILLU_CIE_XW * ( ILLU_CIE_YG - ILLU_CIE_YB ) - \
                            ILLU_CIE_YW * ( ILLU_CIE_XG - ILLU_CIE_XB ) + \
                            ILLU_CIE_XG * ILLU_CIE_YB - ILLU_CIE_XB * ILLU_CIE_YG) )
#define ILLU_CIE_C_GD   ( ( 1.0f / ILLU_CIE_YW ) * \
                          ( ILLU_CIE_XW * ( ILLU_CIE_YB - ILLU_CIE_YR ) - \
                            ILLU_CIE_YW * ( ILLU_CIE_XB - ILLU_CIE_XR ) - \
                            ILLU_CIE_XR * ILLU_CIE_YB + ILLU_CIE_XB * ILLU_CIE_YR ) )
#define ILLU_CIE_C_BD   ( ( 1.0f / ILLU_CIE_YW ) * \
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

// Table of RGB to XYZ (no white balance)
static const float  RGB2XYZ[3][3] =
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
        (1.0f - ILLU_CIE_XR - ILLU_CIE_YR) * ILLU_CIE_C_RD / ILLU_CIE_D,
        (1.0f - ILLU_CIE_XG - ILLU_CIE_YG) * ILLU_CIE_C_GD / ILLU_CIE_D,
        (1.0f - ILLU_CIE_XB - ILLU_CIE_YB) * ILLU_CIE_C_BD / ILLU_CIE_D
	}
};

// Table of XYZ to RGB (no white balance)
static const float  XYZ2RGB[3][3] =
{
	{
	    ( ILLU_CIE_YG - ILLU_CIE_XB - ILLU_CIE_XB * ILLU_CIE_YG + ILLU_CIE_XB * ILLU_CIE_XG )
	    / ILLU_CIE_C_RD,
        ( ILLU_CIE_XB - ILLU_CIE_XG - ILLU_CIE_XB * ILLU_CIE_YG + ILLU_CIE_XG * ILLU_CIE_XB )
        / ILLU_CIE_C_RD,
        ( ILLU_CIE_XG * ILLU_CIE_XB - ILLU_CIE_XB * ILLU_CIE_YG )
        / ILLU_CIE_C_RD
	},
	{
	    ( ILLU_CIE_XB - ILLU_CIE_YR - ILLU_CIE_XB * ILLU_CIE_XR + ILLU_CIE_YR * ILLU_CIE_XB )
	    / ILLU_CIE_C_GD,
        ( ILLU_CIE_XR - ILLU_CIE_XB - ILLU_CIE_XR * ILLU_CIE_XB + ILLU_CIE_XB * ILLU_CIE_YR )
        / ILLU_CIE_C_GD,
        ( ILLU_CIE_XB * ILLU_CIE_YR - ILLU_CIE_XR * ILLU_CIE_XB )
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
	return pow (x, b);		// pow(x, log(bias)/log(0.5)
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
        return ( x * (  6.0 + 0.7662 * x ) / ( 5.9897 + 3.7658 * x ) );
	}

	return log(x + 1);
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

Fl_RGB_Image* fl_imgtk_F2RGB( fl_imgtk_fimg* img )
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
            #pragma omp parallel for
            for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
            {
                buff[ cnt ] = (uchar)( (float)img->pixels[ cnt ] * 255.0f );
            }

            return new Fl_RGB_Image( buff, w, h ,d );
        }
    }

    return NULL;
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
            #pragma omp parallel for
            for( unsigned cnt=0; cnt<(w*h*d); cnt++ )
            {
                float conv = ( img->pixels[ cnt ] > 1.0f ) ? 1.0f : img->pixels[ cnt ];
                buff[ cnt ] = (uchar)( conv * 255.0f + 0.5f );
            }

            return new Fl_RGB_Image( buff, w, h ,d );
        }
    }

    return NULL;
}


// Convert in-place floating point RGB data to Yxy.
bool fl_imgtk_F2YXY( fl_imgtk_fimg* img )
{
    if ( img != NULL )
    {
        #pragma omp parellel for
        for( unsigned cnt=0; cnt<(img->w * img->h); cnt++ )
        {
            float tf[3] = {0};
            float* psrc = &img->pixels[ cnt * img->d ];

            for ( unsigned rpt=0; rpt<3; rpt++ )
            {
                tf[rpt] += RGB2XYZ[rpt][0] * psrc[0];
                tf[rpt] += RGB2XYZ[rpt][1] * psrc[1];
                tf[rpt] += RGB2XYZ[rpt][2] * psrc[2];
            }

            float fW = tf[0] + tf[1] + tf[2];
            float fY = tf[1];

            if ( fY > 0.0f )
            {
                psrc[0] = fY;
                psrc[1] = tf[0] / fW;
                psrc[2] = tf[1] / fW;
            }
            else
            {
                psrc[0] = 0;
                psrc[1] = 0;
                psrc[2] = 0;
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

                float luma790 = ( ( (float)psrc[0] * 0.2126f ) +
                                  ( (float)psrc[1] * 0.7152f ) +
                                  ( (float)psrc[2] * 0.0722f ) );

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
        #pragma omp parellel for
        for( unsigned cnt=0; cnt<(img->w * img->h); cnt++ )
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

            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                tf[rpt] += XYZ2RGB[rpt][0] * psrc[0];
                tf[rpt] += XYZ2RGB[rpt][1] * psrc[1];
                tf[rpt] += XYZ2RGB[rpt][2] * psrc[2];
            }

            memcpy( &psrc, &tf, 3 * sizeof(float) );
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

        maxlumi = 0.0f;
        minlumi = 0.0f;
        worldlumi = 0.0f;

        //#pragma omp parellel for private( dsum )
        for( unsigned cnt=0; cnt<cmax; cnt++ )
        {
            float fY = MAX( 0.0f, img->pixels[ cnt ] );

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

    maxlumi  = -1e20f;
    minlumi  = 1e20f;

    double   lsum  = 0.0;
    double   llsum = 0.0;

    #pragma omp parallel for private( lsum, llsum )
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        float Y = img->pixels[ cnt ];
        maxlumi = ( maxlumi < Y ) ? Y:maxlumi;
        minlumi = ( ( Y > 0.0f ) && ( minlumi < Y ) ) ? minlumi : Y;
        lsum += Y;
        llsum += log( 2.3e-5f + Y );
    }

    avglumi = (float)( lsum / (double)imgsz );
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

    m = ( m > 0.0f ) ? m:(float)( 0.3f + 0.7f * pow(key, 1.4f) );

    float colMax = -1e6f;
    float colMin = +1e6f;

    unsigned imgsz = img->w * img->h;

    if ( ( a == 1.0f ) && ( c == 0.0f ) )
    {
        for( unsigned cnt=0; cnt<imgsz; cnt++ )
        {
            // = pixel luminance
            float pixllumi = Y->pixels[cnt];

            img->pixels[ cnt * img->d ] /= ( img->pixels[ cnt * img->d ]
                                             + pow( f * pixllumi, m ) );

            colMax = ( img->pixels[ cnt * img->d ] > colMax ) ?
                     img->pixels[ cnt * img->d ] : colMax;
            colMin = ( img->pixels[ cnt * img->d ] < colMin ) ?
                     img->pixels[ cnt * img->d ] : colMin;

        }
    }
    else
    {
        float chnlavg[3] = {0.0f};

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
                float apixl    = img->pixels[ cnt * img->d + rpt ];
                // global light adaptation
                float lightloc = c *
                                 + ( 1.0f - c ) * pixllumi;

                // local light adaptation
                float lightglb = c * chnlavg[ rpt ] + ( 1.0f - c ) * lumiavg;

                // interpolated pixel light adaptation
                float lightinp = a * lightloc + ( 1.0f - a ) * lightglb;

                apixl /= apixl + pow( f * lightinp, m );
                img->pixels[ cnt * img->d + rpt ] = apixl;
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
                img->pixels[ cnt * img->d + rpt ] -= colMin;
                img->pixels[ cnt * img->d + rpt ] /= colRange;
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
            double Yw = img->pixels[0] / avgLum;
            Yw *= exposure;
            double interpol = log( 2.0 + fl_imgtk_biasFunction( biasp, Yw / lmax ) * 8.0 );
            double L = fl_imgtk_pade_log( Yw );
            img->pixels[2] = (float)( (L / interpol) / divider);
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
            slope = 4.5f * ( gval - 2.0f ) * 7.5;
        }
        else
        if ( gval <= 1.9f)
        {
            start = 0.018f / ( 2.0f - gval ) * 7.5f;
            slope = 4.5f * ( 2.0f - gval ) * 7.5;
        }

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<(img->w*img->h); cnt++ )
        {
            for( unsigned rpt=0; rpt<img->d; rpt++ )
            {
                float pixel = img->pixels[ cnt * img->d + rpt ];
                pixel = ( pixel <= start ) ? pixel * slope
                        : ( 1.099f * pow( pixel, fgamma) - 0.099f );
                img->pixels[ cnt * img->d + rpt ] = pixel;
            }
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

        if ( fl_imgtk_F2YXY( convimg ) == false )
        {
            fl_imgtk_discard_fimg( convimg );

            return retimg;
        }

        float lumimax = 0.0f;
        float lumimin = 0.0f;
        float lumiavg = 0.0f;

        fl_imgtk_luminancefromYXY( convimg, lumimax, lumimin, lumiavg );
        fl_imgtk_tonemappingdrago( convimg, lumimax, lumiavg, biasparam, exposparam );
        fl_imgtk_YXY2F( convimg );

        if ( gamma != 1.0 )
        {
            fl_imgtk_rec709gammacorrect( convimg, gamma );
        }

        retimg = fl_imgtk_clamp_F2RGB( convimg );

        fl_imgtk_discard_fimg( convimg );
    }

    return retimg;
}
