#include "fl_imgtk.h"
#include "fl_imgtk_minmax.h"

#include <FL/Fl_Image_Surface.H>
#include <FL/fl_draw.H>

#include "fl_smimg.h"

#ifdef USING_OMP
#include <omp.h>
#endif /// of USING_OMP

////////////////////////////////////////////////////////////////////////////////

#define FLOAT_PI        3.141592654

#define fl_imgtk_degree2f( _x_ )      ( ( 360.0 - _x_ ) / 360.0 * 100.0 )
#define fl_imgtk_swap_uc( _a_, _b_ )   uchar t=_a_; _a_=_b_; _b_=t;

////////////////////////////////////////////////////////////////////////////////

const float matrixdata_blur[] =
{
    0.0, 0.2,  0.0,
    0.2, 0.4,  0.2,
    0.0, 0.2,  0.0
};

const float matrixdata_blurmore[] =
{
    0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 2.0, 1.0, 0.0,
    1.0, 2.0, 4.0, 2.0, 1.0,
    0.0, 1.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0
};

const float matrixdata_sharpen[] =
{
    0.0, -1.0,  0.0,
   -1.0,  5.0, -1.0,
    0.0, -1.0,  0.0
};

const float matrixdata_sharpenmore[] =
{
   -1.0, -1.0, -1.0,
   -1.0,  9.0, -1.0,
   -1.0, -1.0, -1.0
};

////////////////////////////////////////////////////////////////////////////////

Fl_RGB_Image* fl_imgtk::fliphorizontal( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        uchar* buff = new uchar[ w * h * d ];

        if ( buff == NULL )
            return NULL;

        memcpy( buff, ptr, w * h * d );

        unsigned hcenter = h/2;
        unsigned cnth = 0;
        unsigned cntw = 0;

        #pragma omp parallel for private(cntw)
        for( cnth=0; cnth<hcenter; cnth++ )
        {
            for( cntw=0; cntw<w; cntw++ )
            {
                unsigned pos1 = ( w * ( h - 1 - cnth ) + cntw ) * d;
                unsigned pos2 = ( w * cnth + cntw ) * d;

                for( unsigned cntd=0; cntd<d; cntd++ )
                {
                    fl_imgtk_swap_uc( buff[ pos1 + cntd ], buff[ pos2 + cntd ] );
                }
            }
        }

        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );

        return newimg;
    }

    return NULL;
}

bool fl_imgtk::fliphorizontal_ex( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return false;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        unsigned hcenter = h/2;
        unsigned cnth = 0;
        unsigned cntw = 0;

        #pragma omp parallel for private(cntw)
        for( cnth=0; cnth<hcenter; cnth++ )
        {
            for( cntw=0; cntw<w; cntw++ )
            {
                unsigned pos1 = ( w * ( h - 1 - cnth ) + cntw ) * d;
                unsigned pos2 = ( w * cnth + cntw ) * d;

                for( unsigned cntd=0; cntd<d; cntd++ )
                {
                    fl_imgtk_swap_uc( ptr[ pos1 + cntd ], ptr[ pos2 + cntd ] );
                }
            }
        }

        img->uncache();

        return true;
    }

    return false;
}

Fl_RGB_Image* fl_imgtk::flipvertical( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) && ( ptr != NULL ) )
    {
        uchar* buff = new uchar[ w * h * d ];

        if ( buff == NULL )
            return NULL;

        memcpy( buff, ptr, w * h * d );

        unsigned wcenter = w/2;
        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<wcenter; cntw++ )
        {
            for( cnth=0; cnth<h; cnth++ )
            {
                unsigned pos1 = ( w * cnth + ( w - cntw ) ) * d;
                unsigned pos2 = ( w * cnth + cntw ) * d;

                for( unsigned cntd=0; cntd<d; cntd++ )
                {
                    fl_imgtk_swap_uc( buff[ pos1 + cntd ], buff[ pos2 + cntd ] );
                }
            }
        }

        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );

        return newimg;
    }

    return NULL;
}

bool fl_imgtk::flipvertical_ex( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return false;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) && ( ptr != NULL ) )
    {
        unsigned wcenter = w/2;
        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<wcenter; cntw++ )
        {
            for( cnth=0; cnth<h; cnth++ )
            {
                unsigned pos1 = ( w * cnth + ( w - cntw ) ) * d;
                unsigned pos2 = ( w * cnth + cntw ) * d;

                for( unsigned cntd=0; cntd<d; cntd++ )
                {
                    fl_imgtk_swap_uc( ptr[ pos1 + cntd ], ptr[ pos2 + cntd ] );
                }
            }
        }

        img->uncache();

        return true;
    }

    return false;
}

Fl_RGB_Image* fl_imgtk::rotate90( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    unsigned src_w = w;
    unsigned src_h = h;

    if ( ( src_w > 0 ) && ( src_h > 0 ) && ( ptr != NULL ) )
    {
        unsigned new_w = src_h;
        unsigned new_h = src_w;

        uchar* buff = new uchar[ new_w * new_h * d ];

        if ( buff == NULL )
            return NULL;

        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private( cnth )
        for( cntw=new_w-1; cntw>0; cntw-- )
        {
            for( cnth=0; cnth<new_h; cnth++ )
            {
                unsigned pos1 = ( new_w * cnth + cntw ) * d;
                unsigned pos2 = ( src_w * ( new_w - cntw - 1 ) + cnth ) * d;

                memcpy( &buff[ pos1 ], &ptr[ pos2 ], d );
            }
        }

        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, new_w, new_h, d );

        return newimg;
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::rotate180( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    unsigned cur_w = w;
    unsigned cur_h = h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) && ( ptr != NULL ) )
    {
        uchar* buff = new uchar[ w * h * d ];

        if ( buff == NULL )
            return NULL;

        memcpy( buff, ptr, w * h * d );

        unsigned imgmax = w*h;
        unsigned cntmax = imgmax / 2;

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<cntmax; cnt++ )
        {
            for( unsigned cntd=0;cntd<d; cntd++)
            {
                fl_imgtk_swap_uc( buff[ cnt * d + cntd ],
                                  buff[ (imgmax - cnt) * d + cntd ] );
            }
        }

        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );

        return newimg;
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::rotate270( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();

    unsigned src_w = w;
    unsigned src_h = h;

    if ( ( src_w > 0 ) && ( src_h > 0 ) )
    {
        unsigned new_w = src_h;
        unsigned new_h = src_w;

        uchar* buff = new uchar[ new_w * new_h * d ];

        if ( buff == NULL )
            return NULL;

        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private( cnth )
        for( cntw=0; cntw<new_w; cntw++ )
        {
            for( cnth=new_h-1; cnth>0; cnth-- )
            {
                unsigned pos1 = ( new_w * cnth + cntw ) * d;
                unsigned pos2 = ( src_w * cntw + new_h - cnth ) * d;

                memcpy( &buff[ pos1 ], &ptr[ pos2 ], d );
            }
        }

        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, new_w, new_h, d );

        return newimg;
    }

    return NULL;
}

float fl_imgtk_min4(float a, float b, float c, float d)
{
   float mn = a;
   if(mn > b) mn = b;
   if(mn > c) mn = c;
   if(mn > d) mn = d;
   return mn;
}

float fl_imgtk_max4(float a, float b, float c, float d)
{
   float mx = a;
   if(mx < b) mx = b;
   if(mx < c) mx = c;
   if(mx < d) mx = d;
   return mx;
}

// rotatefree() Code inspired from
// http://www.codeguru.com/cpp/g-m/gdi/article.php/c3693/Rotate-a-Bitmap-at-Any-Angle-Without-GetPixelSetPixel.htm
// Maybe, it will slowly works !
Fl_RGB_Image* fl_imgtk::rotatefree( Fl_RGB_Image* img, float deg )
{
    if ( ( img == NULL ) || ( img->d() < 3 ) )
        return NULL;

    int img_w = img->w();
    int img_h = img->h();
    int img_d = img->d();

	float CtX = ( (float) img_w ) / 2.0f;
	float CtY = ( (float) img_h ) / 2.0f;

    float fdeg = fl_imgtk_degree2f( deg );

	float cA = (float)cos( fdeg );
	float sA = (float)sin( fdeg );

	float x1 = CtX + (-CtX) * cA - (-CtY) * sA;
	float x2 = CtX + (img_w - CtX) * cA - (-CtY) * sA;
	float x3 = CtX + (img_w - CtX) * cA - (img_h - CtY) * sA;
	float x4 = CtX + (-CtX) * cA - (img_h - CtY) * sA;

	float y1 = CtY + (-CtY) * cA + (-CtX) * sA;
	float y2 = CtY + (img_h - CtY) * cA + (-CtX) * sA;
	float y3 = CtY + (img_h - CtY) * cA + (img_w - CtX) * sA;
	float y4 = CtY + (-CtY) * cA + (img_w - CtX) * sA;

	int OfX = ((int)floor(fl_imgtk_min4(x1, x2, x3, x4)));
	int OfY = ((int)floor(fl_imgtk_min4(y1, y2, y3, y4)));

	int dstW = ((int)ceil(fl_imgtk_max4(x1, x2, x3, x4))) - OfX;
	int dstH = ((int)ceil(fl_imgtk_max4(y1, y2, y3, y4))) - OfY;

    // Now new image !
    uchar* obuff = new uchar[ dstW * dstH * img_d ];

    if ( obuff == NULL )
        return NULL;

    //memset( obuff, 0, dstW * dstH * img_d );

    uchar* psrc = (uchar*)img->data()[0];

    int stepY = 0;
    int stepX = 0;

    // pointer to destination.
    uchar* dst = obuff;

	#pragma omp parellel for private( stepX )
	for ( stepY = 0; stepY<dstH; stepY++ )
	{
		for ( stepX = 0; stepX<dstW; stepX++ )
		{
#if USING_INTERPOLATED_ROTATE_FREE
            float CtX2 = CtX - OfX;
            float CtY2 = CtY - OfY;

            float orgX = ( cA*(stepX-CtX2) + sA*(stepY-CtY2)) + CtX;
            float orgY = (-sA*(stepX-CtX2) + cA*(stepY-CtY2)) + CtY;

			int iorgX  = (int) orgX;
			int iorgY  = (int) orgY;

            float diffX = (orgX - iorgX);
            float diffY = (orgY - iorgY);

            if ((orgX >= 0) && (orgY >= 0) && (orgX < img_w-1) && (orgY < img_h-1))
            {
                uchar* pd = &obuff[ ( stepY * dstW + stepX ) * img_d ];
                uchar* ps = &psrc[ ( iorgX + iorgY * img_w ) * img_d ];

                // Doing interpolated pixel calculation .
                for( unsigned cntd=0; cntd<img_d; cntd++ )
                {
                    float pv[4] = {0.0f};

                    pv[0] = (float)ps[ cntd ];
                    pv[1] = (float)ps[ cntd + img_d ]; /// Right pixel
                    pv[2] = (float)ps[ cntd + ( img_w * img_d ) ]; /// Below pixel
                    pv[3] = (float)ps[ cntd + ( ( img_w + 1 ) * img_d ) ]; /// Right below pixel.

                    float pc = pv[0] * ( 1.0f - diffX ) * ( 1.0f - diffY ) +
                               pv[1] *          diffX   * ( 1.0f - diffY ) +
                               pv[2] * ( 1.0f - diffX ) *          diffY +
                               pv[3] *          diffX   *          diffY;

                    pd[ cntd ] = (uchar)(pc + 0.5f);
                }
            }
#else
            float orgX = CtX + cA * ((float)stepX + OfX - CtX) + sA * ((float)stepY + OfY - CtY);
            float orgY = CtY - sA * ((float)stepX + OfX - CtX) + cA * ((float)stepY + OfY - CtY);
            float nxtY = CtY - sA * ((float)stepX + OfX - CtX) + cA * ((float)stepY + OfY - CtY - 1);

			int iorgX = (int)orgX;
			int iorgY = (int)orgY;
			int inxtY = (int)nxtY;

			if ((iorgX >= 0) && (iorgY >= 0) && (iorgX < img_w) && (iorgY < img_h))
            {
				memcpy( &dst[ stepX * img_d ],
                        &psrc[ ( iorgX + iorgY * img_w ) * img_d ],
                        img_d );
			}
#endif
		}

		dst += dstW * img_d;
	}

    return new Fl_RGB_Image( obuff, dstW, dstH, img->d() );
}

Fl_RGB_Image* fl_imgtk_curve( Fl_RGB_Image* img, const uchar* LUT )
{
    if ( img == NULL )
        return NULL;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();
    unsigned imgsz = w*h;

    if ( imgsz == 0 )
        return NULL;

    uchar* buff = new uchar[ imgsz * d ];

    if ( buff == NULL )
        return NULL;

    memcpy( buff, ptr, imgsz * d );

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        for( unsigned cntd=1; cntd<=d; cntd++ )
        {
            buff[ cnt * d + cntd ] = LUT[ buff[ cnt * d + cntd ] ];
        }
    }

    return new Fl_RGB_Image( buff, w, h, d );
}

bool fl_imgtk_curve_ex( Fl_RGB_Image* img, const uchar* LUT )
{
    if ( img == NULL )
        return false;

    uchar* ptr = (uchar*)img->data()[0];
    unsigned w = img->w();
    unsigned h = img->h();
    unsigned d = img->d();
    unsigned imgsz = w*h;

    if ( imgsz == 0 )
        return false;

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        for( unsigned cntd=1; cntd<=d; cntd++ )
        {
            ptr[ cnt * d + cntd ] = LUT[ ptr[ cnt * d + cntd ] ];
        }
    }

    img->uncache();

    return true;
}

Fl_RGB_Image* fl_imgtk::gamma( Fl_RGB_Image* img, double gamma )
{
    uchar lut[256] = {0};

    double exponent = 1.0f / gamma;
    double val = 255.0 * pow( 255, -exponent );

    for( unsigned cnt=0; cnt<256; cnt++ )
    {
        double col = pow( (double)cnt, exponent ) * val;

        if ( col > 255 )
        {
            col = 255;
        }

        lut[ cnt ] = (uchar)floor( col + 0.5 );
    }

    return fl_imgtk_curve( img, lut );
}

bool fl_imgtk::gamma_ex( Fl_RGB_Image* img, double gamma )
{
    uchar lut[256] = {0};

    double exponent = 1.0f / gamma;
    double val = 255.0 * pow( 255, -exponent );

    for( unsigned cnt=0; cnt<256; cnt++ )
    {
        double col = pow( (double)cnt, exponent ) * val;

        if ( col > 255 )
        {
            col = 255;
        }

        lut[ cnt ] = (uchar)floor( col + 0.5 );
    }

    return fl_imgtk_curve_ex( img, lut );
}

Fl_RGB_Image* fl_imgtk::brightness( Fl_RGB_Image* img, double perc )
{
    uchar lut[256] = {0};

    const double scale = ( 100.0 + perc ) / 100.0;
    double val = 0.0;

    for( unsigned cnt=0; cnt<256; cnt++ )
    {
        val = (double)cnt * scale;
        val = MAX( 0.0, MIN( val, 255.0 ) );
        lut[ cnt ] = (uchar)floor( val + 0.5 );
    }

    return fl_imgtk_curve( img, lut );
}

bool fl_imgtk::brightbess_ex( Fl_RGB_Image* img, double perc )
{
    uchar lut[256] = {0};

    const double scale = ( 100.0 + perc ) / 100.0;
    double val = 0.0;

    for( unsigned cnt=0; cnt<256; cnt++ )
    {
        val = (double)cnt * scale;
        val = MAX( 0.0, MIN( val, 255.0 ) );
        lut[ cnt ] = (uchar)floor( val + 0.5 );
    }

    return fl_imgtk_curve_ex( img, lut );
}

Fl_RGB_Image* fl_imgtk::contrast( Fl_RGB_Image* img, double perc )
{
    uchar lut[256] = {0};

    const double scale = ( 100.0 + perc ) / 100.0;
    double val = 0.0;

    for( unsigned cnt=0; cnt<256; cnt++ )
    {
        val = 128.0 + ( (double)cnt - 128.0 ) * scale;
        val = MAX( 0.0, MIN( val, 255.0 ) );
        lut[ cnt ] = (uchar)floor( val + 0.5 );
    }

    return fl_imgtk_curve( img, lut );
}

bool fl_imgtk::contrast_ex(  Fl_RGB_Image* img, double perc )
{
    uchar lut[256] = {0};

    const double scale = ( 100.0 + perc ) / 100.0;
    double val = 0.0;

    for( unsigned cnt=0; cnt<256; cnt++ )
    {
        val = 128.0 + ( (double)cnt - 128.0 ) * scale;
        val = MAX( 0.0, MIN( val, 255.0 ) );
        lut[ cnt ] = (uchar)floor( val + 0.5 );
    }

    return fl_imgtk_curve_ex( img, lut );
}

Fl_RGB_Image* fl_imgtk::invert( Fl_RGB_Image* img )
{
    if ( img != NULL )
    {
        unsigned img_w = img->w();
        unsigned img_h = img->h();
        unsigned img_d = img->d();

        if ( ( img_w == 0 ) || ( img_h == 0 ) || ( img_d < 3 ) )
            return NULL;

        unsigned buffsz = img_w * img_h;
        uchar* newbuff = new uchar[ buffsz * img_d ];

        if ( newbuff != NULL )
        {
            uchar* refbuff = (uchar*)img->data()[0];

            #pragma omp parallel for
            for( unsigned cnt=0; cnt<buffsz; cnt++ )
            {
                uchar* psrc = &refbuff[ cnt * img_d ];
                uchar* pdst = &newbuff[ cnt * img_d ];

                // alpha channel will skipped to invert.
                for( unsigned rpt=0; rpt<3; rpt++ )
                {
                    pdst[ rpt ] = 0xFF - psrc[ rpt ];
                }
            }

            return new Fl_RGB_Image( newbuff, img_w, img_h, img_d );
        }
    }

    return NULL;
}

bool fl_imgtk::invert_ex( Fl_RGB_Image* img )
{
    if ( img != NULL )
    {
        unsigned img_w = img->w();
        unsigned img_h = img->h();
        unsigned img_d = img->d();

        if ( ( img_w == 0 ) || ( img_h == 0 ) || ( img_d < 3 ) )
            return NULL;

        unsigned buffsz = img_w * img_h;

        uchar* refbuff = (uchar*)img->data()[0];

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<buffsz; cnt++ )
        {
            uchar* psrc = &refbuff[ cnt * img_d ];

            // alpha channel will skipped to invert.
            for( unsigned rpt=0; rpt<3; rpt++ )
            {
                psrc[ rpt ] = 0xFF - psrc[ rpt ];
            }
        }
		
		img->uncache();
		
        return true;
    }
    return false;
}

Fl_RGB_Image* fl_imgtk::filtered( Fl_RGB_Image* img, kfconfig* kfc )
{
    if ( ( img != NULL ) && ( kfc != NULL ) )
    {
        unsigned img_w = img->w();
        unsigned img_h = img->h();
        unsigned img_d = img->d();

        if ( ( img_w == 0 ) || ( img_h == 0 ) || ( img_d < 3 ) )
            return NULL;

        if ( ( kfc->w == 0 ) || ( kfc->h == 0 ) || ( kfc->msz == 0 ) || ( kfc->m == NULL ) )
            return NULL;

        uchar* pixels  = (uchar*)img->data()[0];
        uchar* newbuff = new uchar[ img_w * img_h * img_d ];

        if ( newbuff == NULL )
            return NULL;

        #pragma omp parallel for
        for( unsigned cntx=0; cntx<img_w; cntx++ )
        {
            for( unsigned cnty=0; cnty<img_h; cnty++ )
            {
                double adj[4] = {0.0};

                // -- applying matrix ---
                for( unsigned fcntx=0; fcntx<kfc->w; fcntx++ )
                {
                    for( unsigned fcnty=0; fcnty<kfc->h; fcnty++ )
                    {
                        unsigned posX = ( cntx - kfc->w / 2 + fcntx + img_w )
                                        % img_w;
                        unsigned posY = ( cnty - kfc->h / 2 + fcnty + img_h )
                                        % img_h;

                        unsigned posM = posY * img_w + posX;

                        if ( posM < ( img_w * img_h ) )
                        {
                            for( unsigned cntd=0; cntd<img_d; cntd ++ )
                            {
                                adj[ cntd ] += (double)pixels[ posM * img_d  + cntd ] *
                                               (double)kfc->m[ fcnty * kfc->w + fcntx ];
                            }
                        }
                    }
                }

                for( unsigned cntd=0; cntd<img_d; cntd++ )
                {
                    uchar rpixel = MIN( MAX( kfc->f * adj[ cntd ] + kfc->b, 0 ), 255 );
                    newbuff[ ( cnty * img_w + cntx ) * img_d + cntd ] = rpixel;
                }
            }
        }

        return new Fl_RGB_Image( newbuff, img_w, img_h, img_d );
    }

    return NULL;
}

bool fl_imgtk::filtered_ex( Fl_RGB_Image* img, kfconfig* kfc )
{
    if ( ( img != NULL ) && ( kfc != NULL ) )
    {
        unsigned img_w = img->w();
        unsigned img_h = img->h();
        unsigned img_d = img->d();

        if ( ( img_w == 0 ) || ( img_h == 0 ) || ( img_d < 3 ) )
            return NULL;

        if ( ( kfc->w == 0 ) || ( kfc->h == 0 ) || ( kfc->msz == 0 ) || ( kfc->m == NULL ) )
            return NULL;

        uchar* pixels  = (uchar*)img->data()[0];

        #pragma omp parallel for
        for( unsigned cntx=0; cntx<img_w; cntx++ )
        {
            for( unsigned cnty=0; cnty<img_h; cnty++ )
            {
                double adj[4] = {0.0};

                // -- applying matrix ---
                for( unsigned fcntx=0; fcntx<kfc->w; fcntx++ )
                {
                    for( unsigned fcnty=0; fcnty<kfc->h; fcnty++ )
                    {
                        unsigned posX = ( cntx - kfc->w / 2 + fcntx + img_w )
                                        % img_w;
                        unsigned posY = ( cnty - kfc->h / 2 + fcnty + img_h )
                                        % img_h;

                        unsigned posM = posY * img_w + posX;

                        if ( posM < ( img_w * img_h ) )
                        {
                            for( unsigned cntd=0; cntd<img_d; cntd ++ )
                            {
                                adj[ cntd ] += (double)pixels[ posM * img_d  + cntd ] *
                                               (double)kfc->m[ fcnty * kfc->w + fcntx ];
                            }
                        }
                    }
                }

                for( unsigned cntd=0; cntd<img_d; cntd++ )
                {
                    uchar rpixel = MIN( MAX( kfc->f * adj[ cntd ] + kfc->b, 0 ), 255 );
                    pixels[ ( cnty * img_w + cntx ) * img_d + cntd ] = rpixel;
                }
            }
        }

        return true;
    }

    return false;
}

Fl_RGB_Image* fl_imgtk::rescale( Fl_RGB_Image* img, unsigned w, unsigned h, rescaletype rst )
{
    if ( ( img != NULL ) && ( w > 0 ) && ( h > 0 ) )
    {
        GenericFilter* afilter = NULL;
        Fl_RGB_Image* newimg = NULL;

        switch( (int)rst )
        {
            case (int)NONE:
                afilter = new BoxFilter();
                break;

            case (int)BILINEAR:
                afilter = new BilinearFilter();
                break;

            case (int)BICUBIC:
                afilter = new BicubicFilter();
                break;

            case (int)LANCZOS:
                afilter = new Lanczos3Filter();
                break;

            case (int)BSPLINE:
                afilter = new BSplineFilter();
                break;
        }

        if ( afilter != NULL )
        {
            ResizeEngine* rse = new ResizeEngine( afilter );
            if ( rse != NULL )
            {
                newimg = rse->scale( img, w, h );

                delete rse;
            }

            delete afilter;
        }

        return newimg;
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::draw_widgetimage( Fl_Widget* w )
{
    Fl_RGB_Image* newimg = NULL;

    if ( w != NULL )
    {
        if ( ( w->w() <= 0 ) || ( w->h() <= 0 ) )
        {
            return NULL;
        }

        int prev_visible = w->visible_r();

        if ( prev_visible == 0 )
        {
            w->show();
        }

        Fl_Image_Surface* imgsurf = new Fl_Image_Surface( w->w(), w->h(), 0 );

        if ( imgsurf != NULL )
        {
            imgsurf->set_current();
            imgsurf->draw( w );

            if ( prev_visible == 0 )
            {
                w->hide();
            }

            Fl_Display_Device::display_device()->set_current();

            newimg = imgsurf->image();

            delete imgsurf;
        }
    }

    return newimg;
}

Fl_RGB_Image* fl_imgtk::drawblurred_widgetimage( Fl_Widget* w, unsigned factor )
{
    Fl_RGB_Image* blurredimg = NULL;

    if ( w != NULL )
    {
        if ( ( w->w() <= 0 ) || ( w->h() <= 0 ) )
        {
            return NULL;
        }

        int prev_visible = w->visible_r();

        if ( prev_visible == 0 )
        {
            w->show();
        }

        Fl_Image_Surface* imgsurf = new Fl_Image_Surface( w->w(), w->h(), 0 );

        if ( imgsurf != NULL )
        {
            imgsurf->set_current();
            imgsurf->draw( w );

            if ( prev_visible == 0 )
            {
                w->hide();
            }

            Fl_Display_Device::display_device()->set_current();

            // Calc scaling factor.
            if ( factor ==  0 )
                factor = 1;

            unsigned scd_w = w->w() / factor;
            unsigned scd_h = w->h() / factor;

            if ( scd_w == 0 )
                scd_w = 10;

            if ( scd_h == 0 )
                scd_h = 10;

            Fl_RGB_Image* widgetimg = imgsurf->image();

            BilinearFilter* blfilter = new BilinearFilter();
            if ( blfilter != NULL )
            {
                ResizeEngine* redown = new ResizeEngine( blfilter );
                if ( redown != NULL )
                {
                    Fl_RGB_Image* sdimg = redown->scale( widgetimg, scd_w, scd_h );
                    if ( sdimg != NULL )
                    {
                        BSplineFilter* bsfilter = new BSplineFilter();
                        if ( bsfilter != NULL )
                        {
                            ResizeEngine* reup = new ResizeEngine( bsfilter );
                            if (  reup != NULL )
                            {
                                blurredimg = reup->scale( sdimg, w->w(), w->h() );

                                delete reup;
                            }
                            delete bsfilter;
                        }
                        discard_user_rgb_image( sdimg );
                    }
                    delete redown;
                }
                delete blfilter;
            }

            delete widgetimg;
        }
    }

    return blurredimg;
}

Fl_RGB_Image* fl_imgtk::blurredimage( Fl_RGB_Image* src, unsigned factor )
{
    Fl_RGB_Image* newimg = NULL;

    if ( src != NULL )
    {
        // Calc scaling factor.
        if ( factor ==  0 )
            factor = 1;

        unsigned scd_w = src->w() / factor;
        unsigned scd_h = src->h() / factor;

        if ( scd_w == 0 )
            scd_w = 10;

        if ( scd_h == 0 )
            scd_h = 10;

        BilinearFilter* blfilter = new BilinearFilter();
        if ( blfilter != NULL )
        {
            ResizeEngine* redown = new ResizeEngine( blfilter );
            if ( redown != NULL )
            {
                Fl_RGB_Image* sdimg = redown->scale( src, scd_w, scd_h );
                if ( sdimg != NULL )
                {
                    BSplineFilter* bsfilter = new BSplineFilter();
                    if ( bsfilter != NULL )
                    {
                        ResizeEngine* reup = new ResizeEngine( bsfilter );
                        if (  reup != NULL )
                        {
                            newimg = reup->scale( sdimg, src->w(), src->h() );

                            delete reup;
                        }
                        delete bsfilter;
                    }
                    discard_user_rgb_image( sdimg );
                }
                delete redown;
            }
            delete blfilter;
        }
    }

    return newimg;
}

bool fl_imgtk::blurredimage_ex( Fl_RGB_Image* src, unsigned factor )
{
    if ( src != NULL )
    {
        // Calc scaling factor.
        if ( factor ==  0 )
            factor = 1;

        unsigned scd_w = src->w() / factor;
        unsigned scd_h = src->h() / factor;

        if ( scd_w == 0 )
            scd_w = 10;

        if ( scd_h == 0 )
            scd_h = 10;

        bool retb = false;

        BilinearFilter* blfilter = new BilinearFilter();
        if ( blfilter != NULL )
        {
            ResizeEngine* redown = new ResizeEngine( blfilter );
            if ( redown != NULL )
            {
                Fl_RGB_Image* sdimg = redown->scale( src, scd_w, scd_h );
                if ( sdimg != NULL )
                {
                    BSplineFilter* bsfilter = new BSplineFilter();
                    if ( bsfilter != NULL )
                    {
                        ResizeEngine* reup = new ResizeEngine( bsfilter );
                        if (  reup != NULL )
                        {
                            Fl_RGB_Image* newimg = reup->scale( sdimg, src->w(), src->h() );

                            if ( newimg != NULL )
                            {
                                uchar* pdst = (uchar*)src->data()[0];
                                uchar* psrc = (uchar*)newimg->data()[0];

                                memcpy( pdst, psrc, src->w() * src->h() * src->d() );
                                src->uncache();

                                discard_user_rgb_image( newimg );

                                retb = true;
                            }

                            delete reup;
                        }
                        delete bsfilter;
                    }
                    discard_user_rgb_image( sdimg );
                }
                delete redown;
            }
            delete blfilter;
        }

        return retb;
    }

    return false;
}

Fl_RGB_Image* fl_imgtk::crop( Fl_RGB_Image* src, unsigned sx, unsigned sy, unsigned w, unsigned h )
{
    if ( src != NULL )
    {
        unsigned rsx = sx;
        unsigned rsy = sy;
        unsigned rw  = w;
        unsigned rh  = h;
        unsigned sd  = src->d();

        if ( ( rsx > src->w() ) || ( rsy > src->h() ) )
            return NULL;

        if ( src->w() > ( rw + rsx ) )
        {
            rw = src->w() - rsx;
        }

        if ( src->h() > ( rh + rsy ) )
        {
            rh = src->h() - rsy;
        }

        uchar* rbuff = (uchar*)src->data()[0];
        uchar* obuff = new uchar[ rw * rh * sd ];

        if ( ( rbuff != NULL ) && ( obuff != NULL ) )
        {
            unsigned srcw = src->w();
            unsigned srch = src->h();
            unsigned dsth = srch - rsy;
            unsigned cnty;
            unsigned cntx;

            #pragma omp parellel for private( cntx )
            for( cnty=rsy; cnty<srch; cnty++ )
            {
                for( cntx=rsx; cntx<srcw; cntx ++ )
                {
                    uchar* rptr = &rbuff[ ( cnty * srcw + cntx ) * sd ];
                    uchar* wptr = &obuff[ ( (cnty - rsy) * rw + ( cntx - rsx) ) * sd ];

                    memcpy( wptr, rptr, sd );
                }
            }

            return new Fl_RGB_Image( obuff, rw, rh , sd );
        }
    }

    return NULL;
}

void fl_imgtk_putimgonbuffer( uchar* buff, unsigned bw, unsigned bh, unsigned bd,
                              Fl_RGB_Image* img, int px, int py, float alpha )
{
    if ( ( buff != NULL ) && ( img != NULL ) )
    {
        if ( img->d() < 3 )
            return;

        uchar* rbuff = (uchar*)img->data()[0];

        int minx = px;
        int miny = py;
        int imgx = 0;
        int imgy = 0;

        int maxw = MIN( img->w() + minx, bw );
        int maxh = MIN( img->h() + miny, bh );

        if ( px < 0 )
        {
            minx = 0;
            maxw = MIN( img->w() + px, bw );
            imgx = abs( px );
            px   = 0;
        }

        if ( py < 0 )
        {
            miny = 0;
            maxh = MIN ( img->h() + py, bh );
            imgy = abs( py );
            py   = 0;
        }

        int cntx = 0;
        int cnty = 0;
        int imgw = img->w();
        int imgd = img->d();

        #pragma omp parellel for private( cnty )
        for( cnty=miny; cnty<maxh; cnty++ )
        {
            for( cntx=minx; cntx<maxw; cntx++ )
            {
                unsigned ipx = imgx + (cntx - minx);
                unsigned ipy = imgy + (cnty - miny);

                uchar* rptr = &rbuff[ ( ipy * imgw + ipx ) * imgd ];
                uchar* wptr = &buff[ ( ( cnty * bw ) + cntx ) * bd ];

                uchar alp = 255;

                if ( imgd == 4 )
                {
                    alp = rptr[3];
                }

                float falp = ( (float)alp / 255.0f ) * alpha;

                for( unsigned rpt=0; rpt<3; rpt++ )
                {
                    float fp = (float)wptr[ rpt ];

                    // Check has alpha ...
                    if ( falp > 0.0f )
                    {
                        // Scale down previous pixel.
                        fp *= ( 1.0f - falp );
                        // And now add new alpha pixel.
                        fp += (float)rptr[ rpt ] * falp;

                        // Cutoff maximum.
                        fp = MIN( 255.0f, fp );
                    }

                    wptr[ rpt ] = (uchar)fp;

                }

                if ( bd == 4 )
                {
                    wptr[3] = 0xFF;
                }
            }
        }
    }
}

void fl_imgtk_subimgonbuffer( uchar* buff, unsigned bw, unsigned bh, unsigned bd,
                              Fl_RGB_Image* img, int px, int py, float alpha )
{
    if ( ( buff != NULL ) && ( img != NULL ) )
    {
        if ( img->d() < 3 )
            return;

        uchar* rbuff = (uchar*)img->data()[0];

        int minx = px;
        int miny = py;
        int imgx = 0;
        int imgy = 0;

        int maxw = MIN( img->w() + minx, bw );
        int maxh = MIN( img->h() + miny, bh );

        if ( px < 0 )
        {
            minx = 0;
            maxw = MIN( img->w() + px, bw );
            imgx = abs( px );
            px   = 0;
        }

        if ( py < 0 )
        {
            miny = 0;
            maxh = MIN ( img->h() + py, bh );
            imgy = abs( py );
            py   = 0;
        }

        int cntx = 0;
        int cnty = 0;
        int imgw = img->w();
        int imgd = img->d();

        #pragma omp parellel for private( cnty )
        for( cnty=miny; cnty<maxh; cnty++ )
        {
            for( cntx=minx; cntx<maxw; cntx++ )
            {
                unsigned ipx = imgx + (cntx - minx);
                unsigned ipy = imgy + (cnty - miny);

                uchar* rptr = &rbuff[ ( ipy * imgw + ipx ) * imgd ];
                uchar* wptr = &buff[ ( ( cnty * bw ) + cntx ) * bd ];

                uchar alp = 255;

                if ( imgd == 4 )
                {
                    alp = rptr[3];
                }

                float falp = ( (float)alp / 255.0f ) * alpha;

                for( unsigned rpt=0; rpt<3; rpt++ )
                {
                    float fp = (float)wptr[ rpt ];

                    // Check has alpha ...
                    if ( falp > 0.0f )
                    {
                        // Scale down previous pixel.
                        fp *= ( 1.0f - falp );
                        // And now subtract new alpha pixel.
                        fp -= (float)rptr[ rpt ] * falp;

                        // Cutoff maximum.
                        fp = MAX( 0.0f, fp );
                    }

                    wptr[ rpt ] = (uchar)fp;

                }

                if ( bd == 4 )
                {
                    wptr[3] = 0xFF;
                }
            }
        }
    }
}


Fl_RGB_Image* fl_imgtk::merge( Fl_RGB_Image* src1, Fl_RGB_Image* src2, mergeconfig* cfg )
{
    Fl_RGB_Image* newimg = NULL;

    if ( ( src1 != NULL ) && ( src2 != NULL ) )
    {
        float fratios[2] = {0};

        unsigned maxsz_w = src1->w();
        unsigned maxsz_h = src1->h();
        int      img1px  = 0;
        int      img1py  = 0;
        int      img2px  = 0;
        int      img2py  = 0;

        // Recognize working conditions  ....
        if ( cfg != NULL )
        {
            fratios[0] = MIN( 1.0f, MAX( 0.0f, cfg->src1ratio ) );
            fratios[1] = MIN( 1.0f, MAX( 0.0f, cfg->src2ratio ) );

            if ( cfg->autoexpand == true )
            {
                int minx = MIN( cfg->src1putx, cfg->src2putx );
                int miny = MIN( cfg->src1puty, cfg->src2puty );
                int maxw = MAX( src1->w(), src2->w() );
                int maxh = MAX( src1->h(), src2->h() );

                maxsz_w = maxw - minx;
                maxsz_h = maxh - miny;
            }

            img1px = cfg->src1putx;
            img1py = cfg->src1puty;
            img2px = cfg->src2putx;
            img2py = cfg->src2puty;
        }
        else
        {
            fratios[0] = 1.0f;
            fratios[1] = 1.0f;
        }

        uchar* obuff = new uchar[ maxsz_w * maxsz_h * 4 ];

        if ( obuff == NULL )
            return  NULL;

        memset( obuff, 0, maxsz_w * maxsz_h * 4 );

        fl_imgtk_putimgonbuffer( obuff, maxsz_w, maxsz_h, 4,
                                 src1, img1px, img1py, fratios[0] );
        fl_imgtk_putimgonbuffer( obuff, maxsz_w, maxsz_h, 4,
                                 src2, img2px, img2py, fratios[1] );

        newimg = new Fl_RGB_Image( obuff, maxsz_w, maxsz_h, 4 );
    }

    return newimg;
}

Fl_RGB_Image* fl_imgtk::subtract( Fl_RGB_Image* src1, Fl_RGB_Image* src2, int px, int py, float sr )
{
    Fl_RGB_Image* newimg = NULL;

    if ( ( src1 != NULL ) && ( src2 != NULL ) )
    {
        unsigned maxsz_w = src1->w();
        unsigned maxsz_h = src1->h();

        // Recognize working conditions  ....

        uchar* obuff = new uchar[ maxsz_w * maxsz_h * 4 ];

        if ( obuff == NULL )
            return  NULL;

        memset( obuff, 0, maxsz_w * maxsz_h * 4 );

        fl_imgtk_putimgonbuffer( obuff, maxsz_w, maxsz_h, 4, src1, 0, 0, 1.0f );
        fl_imgtk_subimgonbuffer( obuff, maxsz_w, maxsz_h, 4, src2, px, py, sr );

        newimg = new Fl_RGB_Image( obuff, maxsz_w, maxsz_h, 4 );
    }

    return newimg;
}

bool fl_imgtk::sbutract_ex( Fl_RGB_Image* src1, Fl_RGB_Image* src2, int px, int py, float sr )
{
    if ( ( src1 != NULL ) && ( src2 != NULL ) )
    {
        uchar* obuff = (uchar*)src1->data()[0];

        fl_imgtk_subimgonbuffer( obuff, src1->w(), src2->w(), src1->d(),
                                 src2, px, py, sr );

        src1->uncache();

        return true;
    }

    return false;
}


unsigned fl_imgtk::makealphamap( uchar* &amap, Fl_RGB_Image* src, float val )
{
    if ( src == NULL )
        return 0;

    if ( ( src->w() == 0 ) || ( src->h() == 0 ) || ( src->d() < 3 ) )
        return 0;

    val = MAX( 1.0f, MIN( 0.0f, val ) );

    unsigned imgsz = src->w() * src->h();

    uchar* refbuff = (uchar*)src->data()[0];
    uchar* newbuff = new uchar[ imgsz ];

    if ( newbuff != NULL )
    {
        if ( src->d() == 4 )
        {
            #pragma omp parellel for
            for( unsigned cnt=0; cnt<imgsz; cnt++ )
            {
                uchar fillv = uchar( (float)refbuff[cnt * 4] * val );

                newbuff[ cnt ] = fillv;
            }
        }
        else
        {
            // Generate alphamap from RGB average.

            #pragma omp parellel for
            for( unsigned cnt=0; cnt<imgsz; cnt++ )
            {
                unsigned avrp = refbuff[cnt*3+0] + refbuff[cnt*3+1] + refbuff[cnt*3+2];
                avrp /= 3;

                newbuff[ cnt ] = uchar( (float)avrp * val );
            }
        }

        amap = newbuff;
        return imgsz;
    }

    return 0;
}

unsigned fl_imgtk::makealphamap( uchar* &amap, unsigned w, unsigned h, uchar val )
{
    unsigned imgsz = w * h;

    val = MAX( 1.0f, MIN( 0.0f, val ) );

    uchar* newbuff = new uchar[ imgsz ];

    if ( newbuff != NULL )
    {
        uchar fillv = (uchar)( 255.0f * val );

        memset( newbuff, fillv, imgsz );

        amap = newbuff;
        return imgsz;
    }

    return 0;
}

Fl_RGB_Image* fl_imgtk::applyalpha( Fl_RGB_Image* src, float val )
{
    Fl_RGB_Image* newimg = NULL;

    if ( src != NULL )
    {
        if ( src->d() < 3 )
            return NULL;

        unsigned img_w  = src->w();
        unsigned img_h  = src->h();
        unsigned img_sz = img_w * img_h;

        uchar* sbuff = (uchar*)src->data()[0];
        uchar* obuff = new uchar[ img_sz * 4 ];

        if ( obuff == NULL )
            return NULL;

        val = MIN( 1.0f, MAX( 0.0f, val ) );

        // Copy pixels from source image.
        if ( src->d() == 4 )
        {
            memcpy( obuff, sbuff, img_sz * 4 );
        }
        else
        {
            #pragma omp parellel for
            for( unsigned cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* sbpos = &sbuff[ cnt * 3 ];
                uchar* obpos = &obuff[ cnt * 4 ];

                memcpy( obpos , sbpos, 3 );
                obpos[3] = 0xFF;
            }
        }

        #pragma omp parellel for
        for( unsigned cnt=0; cnt<img_sz; cnt++ )
        {
            uchar* obp = &obuff[ cnt * 4 + 3 ];

            *obp = (uchar)( (float)*obp * val );
        }
		
		newimg = new Fl_RGB_Image( obuff, img_w, img_h, 4 );
    }

    return newimg;
}

Fl_RGB_Image* fl_imgtk::applyalpha( Fl_RGB_Image* src, uchar* alphamap, unsigned amsz )
{
    Fl_RGB_Image* newimg = NULL;

    if ( src != NULL )
    {
        if ( src->d() < 3 )
            return NULL;

        unsigned img_w  = src->w();
        unsigned img_h  = src->h();
        unsigned img_sz = img_w * img_h;

        uchar* sbuff = (uchar*)src->data()[0];
        uchar* obuff = new uchar[ img_sz * 4 ];

        if ( obuff == NULL )
            return NULL;

        // Copy pixels from source image.
        if ( src->d() == 4 )
        {
            memcpy( obuff, sbuff, img_sz * 4 );
        }
        else
        {
            #pragma omp parellel for
            for( unsigned cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* sbpos = &sbuff[ cnt * 3 ];
                uchar* obpos = &obuff[ cnt * 4 ];

                memcpy( obpos , sbpos, 3 );
            }
        }

        // Apply alpha map.
        if ( ( alphamap != NULL ) && ( amsz == img_sz ) )
        {
            #pragma omp parellel for
            for( unsigned cnt=0; cnt<img_sz; cnt++ )
            {
                obuff[ cnt * 4 + 3 ] = alphamap[ cnt ];
            }
        }

        newimg = new Fl_RGB_Image( obuff, img_w, img_h, 4 );
    }

    return newimg;
}

bool fl_imgtk::applyalpha_ex( Fl_RGB_Image* src, float val )
{
    if ( src != NULL )
    {
        if ( src->d() < 3 )
            return NULL;

        unsigned img_w  = src->w();
        unsigned img_h  = src->h();
        unsigned img_sz = img_w * img_h;

        uchar* ptr = (uchar*)src->data()[0];

        val = MAX( 1.0f, MIN( 0.0f, val ) );

        if ( src->d() == 4 )
        {
            #pragma omp parellel for
            for( unsigned cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* obp = &ptr[ cnt * 4 + 3 ];

                *obp = (uchar)( (float)*obp * val );
            }
        }
        else
        {
            #pragma omp parellel for
            for( unsigned cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* obp = &ptr[ cnt * 3 ];

                for( unsigned rpt=0; rpt<3; rpt++ )
                {
                    obp[rpt] = (uchar)( (float)obp[rpt] * val );
                }
            }
        }
    }

    return false;
}

bool fl_imgtk::drawonimage( Fl_RGB_Image* bgimg, Fl_RGB_Image* img, int x, int y, float alpha )
{
    if ( ( bgimg != NULL ) && ( img != NULL ) )
    {
        uchar* bgbuff = (uchar*)bgimg->data()[0];

        fl_imgtk_putimgonbuffer( bgbuff, bgimg->w(), bgimg->h(), bgimg->d(),
                                 img, x, y, alpha );

        bgimg->uncache();

        return true;
    }

    return false;
}

Fl_RGB_Image* fl_imgtk::makegradation_h( unsigned w, unsigned h, ulong col1, ulong col2, bool dither )
{
    uchar* buffer = new unsigned char[ w * h * 3 ];
    if ( buffer != NULL )
    {
        float downscale_f = 2559.0f / ( float( h ) * 10.0f );
        float col1g     = 0.0f;
        float col2g     = 0.0f;
        float randfv    = 0.0f;
        uchar fill_r    = 0;
        uchar fill_g    = 0;
        uchar fill_b    = 0;
        uchar ref_c1r   = ( col1 & 0xFF000000 ) >> 24;
        uchar ref_c1g   = ( col1 & 0x00FF0000 ) >> 16;
        uchar ref_c1b   = ( col1 & 0x0000FF00 ) >> 8;
        uchar ref_c2r   = ( col2 & 0xFF000000 ) >> 24;
        uchar ref_c2g   = ( col2 & 0x00FF0000 ) >> 16;
        uchar ref_c2b   = ( col2 & 0x0000FF00 ) >> 8;

        #pragma omp parellel for
        for ( int cy=0; cy<h; cy++ )
        {
            uchar* sbuf = &buffer[ cy * w * 3 ];

            col1g = ( ( float( h*10 ) - float( cy*10 )  ) * downscale_f ) / 2559.0f;
            col2g = 1.0f - col1g;

            fill_r = ( (float)ref_c1r * col1g ) + ( (float)ref_c2r * col2g );
            fill_g = ( (float)ref_c1g * col1g ) + ( (float)ref_c2g * col2g );
            fill_b = ( (float)ref_c1b * col1g ) + ( (float)ref_c2b * col2g );

            for ( int cx=0; cx<w; cx++ )
            {
                if ( dither == true )
                {
                    if ( h > 255 )
                    {
                        // A simple trick to ...
                        // Make dithered color !
                        randfv = float( rand() % 256 ) / 5120.0f;
                    }

                    sbuf[0] = fill_r - uchar( ref_c2r * col2g * randfv );
                    sbuf[1] = fill_g - uchar( ref_c2g * col2g * randfv );
                    sbuf[2] = fill_b - uchar( ref_c2b * col2g * randfv );
                }
                else
                {
                    sbuf[0] = fill_r;
                    sbuf[1] = fill_g;
                    sbuf[2] = fill_b;
                }
                sbuf += 3;
            }
        }

        return new Fl_RGB_Image( buffer, w, h, 3 );
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::makegradation_v( unsigned w, unsigned h, ulong col1, ulong col2, bool dither )
{
    uchar* buffer = new unsigned char[ w * h * 3 ];
    if ( buffer != NULL )
    {
        float downscale_f = 2559.0f / ( float( h ) * 10.0f );
        float col1g     = 0.0f;
        float col2g     = 0.0f;
        float randfv    = 0.0f;
        uchar fill_r    = 0;
        uchar fill_g    = 0;
        uchar fill_b    = 0;
        uchar ref_c1r   = ( col1 & 0xFF000000 ) >> 24;
        uchar ref_c1g   = ( col1 & 0x00FF0000 ) >> 16;
        uchar ref_c1b   = ( col1 & 0x0000FF00 ) >> 8;
        uchar ref_c2r   = ( col2 & 0xFF000000 ) >> 24;
        uchar ref_c2g   = ( col2 & 0x00FF0000 ) >> 16;
        uchar ref_c2b   = ( col2 & 0x0000FF00 ) >> 8;

        #pragma omp parellel for
        for ( int cx=0; cx<w; cx++ )
        {
            for ( int cy=0; cy<h; cy++ )
            {
                uchar* sbuf = &buffer[ cy * w * 3 ];

                col1g = ( ( float( h*10 ) - float( cy*10 )  ) * downscale_f ) / 2559.0f;
                col2g = 1.0f - col1g;

                fill_r = ( (float)ref_c1r * col1g ) + ( (float)ref_c2r * col2g );
                fill_g = ( (float)ref_c1g * col1g ) + ( (float)ref_c2g * col2g );
                fill_b = ( (float)ref_c1b * col1g ) + ( (float)ref_c2b * col2g );

                if ( dither == true )
                {
                    if ( h > 255 )
                    {
                        // A simple trick to ...
                        // Make dithered color !
                        randfv = float( rand() % 256 ) / 5120.0f;
                    }

                    sbuf[0] = fill_r - uchar( ref_c2r * col2g * randfv );
                    sbuf[1] = fill_g - uchar( ref_c2g * col2g * randfv );
                    sbuf[2] = fill_b - uchar( ref_c2b * col2g * randfv );
                }
                else
                {
                    sbuf[0] = fill_r;
                    sbuf[1] = fill_g;
                    sbuf[2] = fill_b;
                }
                sbuf += 3;
            }
        }

        return new Fl_RGB_Image( buffer, w, h, 3 );
    }

    return NULL;
}


void fl_imgtk_t_set_kfconfig( fl_imgtk::kfconfig* kfc, uchar w, uchar h, float f, float b, const float* m )
{
    if ( ( kfc != NULL ) && ( m != NULL ) )
    {
        kfc->w = w;
        kfc->h = h;
        kfc->f = f;
        kfc->b = b;
        kfc->msz = w*h;
        kfc->m = new float[ kfc->msz ];
        if ( kfc->m != NULL )
        {
            memcpy( kfc->m, m, kfc->msz * sizeof( float ) );
        }
        else
        {
            kfc->msz = 0;
        }
    }
}

fl_imgtk::kfconfig* fl_imgtk::new_kfconfig( const char* preset )
{
    fl_imgtk::kfconfig* newcfg = new fl_imgtk::kfconfig;

    if ( newcfg != NULL )
    {
        memset( newcfg, 0, sizeof( fl_imgtk::kfconfig ) );

        // check preset conf.
        if ( preset != NULL )
        {
            if ( strcmp( preset, "blur" ) == 0 )
            {
                fl_imgtk_t_set_kfconfig( newcfg,
                                         3, 3, 0.8f, 0.0f,
                                         matrixdata_blur );
            }
            else
            if ( strcmp( preset, "blurmore" ) == 0 )
            {
                fl_imgtk_t_set_kfconfig( newcfg,
                                         5, 5, 1.0f/20.0f, 0.0f,
                                         matrixdata_blurmore );
            }
            else
            if ( strcmp( preset, "sharpen" ) == 0 )
            {
                fl_imgtk_t_set_kfconfig( newcfg,
                                         3, 3, 1.0f, 0.0f,
                                         matrixdata_sharpen );
            }
            else
            if ( strcmp( preset, "sharpenmore" ) == 0 )
            {
                fl_imgtk_t_set_kfconfig( newcfg,
                                         3, 3, 1.0f, 0.0f,
                                         matrixdata_sharpenmore );
            }
        }
    }

    return newcfg;
}

void fl_imgtk::discard_user_rgb_image( Fl_RGB_Image* &img )
{
    if( img != NULL )
    {
        if ( ( img->array != NULL ) && ( img->alloc_array == 0 ) )
        {
            delete[] img->array;
        }

        delete img;
        img = NULL;
    }
}

void fl_imgtk::discard_kfconfig( kfconfig* &kfc )
{
    if ( kfc != NULL )
    {
        if ( ( kfc->msz > 0 ) && ( kfc->m != NULL ) )
        {
            delete[] kfc->m;
        }

        delete kfc;
        kfc = NULL;
    }
}
