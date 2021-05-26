#ifdef _MSC_VER
    // Holy mother F-hawk M$VC ...
    #pragma warning(disable : 4018)
    #pragma warning(disable : 4068)
    #pragma warning(disable : 4244)
    #pragma warning(disable : 4996)
#endif

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint> /// prevent error on Linux.

#include "fl_imgtk.h"
#include "fl_imgtk_minmax.h"
#include "fl_smimg.h"

#include <FL/Fl_Widget.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Image_Surface.H>
#include <FL/fl_draw.H>

using namespace std;

#include "fl_imgtk_cfg.h"

////////////////////////////////////////////////////////////////////////////////

#define FLOAT_PI            3.141592654f
#define FLOAT_PI2X          6.28318530718f

#define FLIMGTK_BI_RGB       0  /// No compression - straight BGR data
#define FLIMGTK_BI_RLE8      1  /// 8-bit run-length compression
#define FLIMGTK_BI_RLE4      2  /// 4-bit run-length compression
#define FLIMGTK_BI_BITFIELDS 3  /// RGB bitmap with RGB masks

#define __MIN(a,b) (((a)<(b))?(a):(b))
#define __MAX(a,b) (((a)>(b))?(a):(b))

#define fl_imgtk_degree2f( _x_ )            ( ( _x_ / 360.f ) * FLOAT_PI2X )
#define fl_imgtk_swap_uc( _a_, _b_ )        uchar _t_=_a_; _a_=_b_; _b_=t

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

inline void fl_imgtk_swap_mem( uchar* a, uchar* b, size_t c )  
{ 
    if ( c > 0 )
    {
        uchar* t = new uchar[c];
        if ( t != NULL )
        {
            memcpy( t, a, c );
            memcpy( a, b, c );
            memcpy( b, t, c );

            delete[] t;
        }
    }
}

Fl_RGB_Image* fl_imgtk::makeanempty( unsigned w, unsigned h, unsigned d, ulong color )
{
    if ( ( w == 0 ) || ( h == 0 ) )
        return NULL;

    uchar ref_r   = ( color & 0xFF000000 ) >> 24;
    uchar ref_g   = ( color & 0x00FF0000 ) >> 16;
    uchar ref_b   = ( color & 0x0000FF00 ) >> 8;
    uchar ref_a   = ( color & 0x000000FF );

    if ( d == 1 )
    {
        float     colaf  = ((float)( ref_r + ref_g + ref_b ) / 3.f) * ( ref_a / 255.f );
        OMPSIZE_T resz   = w * h;
        uchar*    pdata  = new uchar[ resz ];
        uchar     ref_y  = (uchar)colaf;
        
        if ( pdata != NULL )
        {
            memset( pdata, ref_y, resz );
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            Fl_RGB_Image* newimg = new Fl_RGB_Image( pdata, w, h, d );
            if ( newimg != NULL )
            {
                newimg->alloc_array = 1;
                return newimg;
            }
#else
            return new Fl_RGB_Image( pdata, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        }
    }
    else
    if ( d == 2 )
    {
        OMPSIZE_T resz   = w * h;
        OMPSIZE_T datasz = resz * d;
        uchar*    pdata  = new uchar[ datasz ];
        
        ulong colaf   = (float)( ref_r + ref_g + ref_b ) / 3.f;
        uchar ref_y   = (uchar)(colaf & 0x000000FF);

        uchar carray[2] = { ref_y, ref_a };
        
        if ( pdata != NULL )
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<resz; cnt++ )
            {
                memcpy( &pdata[ cnt * d ], &carray[0], d );
            }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            Fl_RGB_Image* newimg = new Fl_RGB_Image( pdata, w, h, d );
            if ( newimg != NULL )
            {
                newimg->alloc_array = 1;
                return newimg;
            }
#else
            return new Fl_RGB_Image( pdata, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        }
    }
    else
    if ( d >= 3 )
    {
        OMPSIZE_T resz   = w * h;
        OMPSIZE_T datasz = resz * d;
        uchar*    pdata  = new uchar[ datasz ];
        
        uchar carray[4] = { ref_r, ref_g, ref_b, ref_a };
        
        if ( pdata != NULL )
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<resz; cnt++ )
            {
                memcpy( &pdata[ cnt * d ], &carray[0], d );
            }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            Fl_RGB_Image* newimg = new Fl_RGB_Image( pdata, w, h, d );
            if ( newimg != NULL )
            {
                newimg->alloc_array = 1;
                return newimg;
            }
#else
            return new Fl_RGB_Image( pdata, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        }
    }
    
    return NULL;
}

inline uint16_t flimgtk_memread_word( const char* buffer, unsigned* que = NULL )
{
    uint16_t retus = 0;
    memcpy( &retus, buffer, 2 );
    
    if ( que != NULL )
    {
        *que += 2;
    }
    
    return retus;
}

inline uint32_t flimgtk_memread_dword( const char* buffer, unsigned* que = NULL )
{
    uint32_t retui = 0;
    memcpy( &retui, buffer, 4 );
    
    if ( que != NULL )
    {
        *que += 4;
    }
    
    return retui;
}

inline int flimgtk_memread_int( const char* buffer, unsigned* que = NULL )
{
    int reti = 0;
    memcpy( &reti, buffer, 4 );
    
    if ( que != NULL )
    {
        *que += 4;
    }
    
    return reti;
}

// Reading BMP code from Fl_BMP_Image.cxx
Fl_RGB_Image* fl_imgtk::createBMPmemory( const char* buffer, unsigned buffersz )
{
    if ( ( buffer != NULL ) && ( buffersz > 32 ) )
    {
        unsigned buffque = 0;
        
        // Check 'BM'.
        if ( ( buffer[buffque++] != 'B' ) || ( buffer[buffque++] != 'M' ) )
            return NULL;
        
        // Skips uselesses.
        buffque += 4;           /// skip size.
        buffque += 2;           /// skip reserved something.
        buffque += 2;
        
        int info_size;          /// Size of info header
        uint16_t depth;         /// Depth of image (bits)
        int bDepth = 3;         /// Depth of image (bytes)
        uint16_t compression;   /// Type of compression
        uint32_t colors_used;   /// Number of colors used
        int x, y;               /// Looping vars
        int32_t color = 0;      /// Color of RLE pixel
        int repcount;           /// Number of times to repeat
        int temp;               /// Temp. Color or Index
        int align;              /// Alignment bytes
        uint32_t dataSize;      /// number of bytes in image data set
        int row_order=-1;       /// 1 = normal;  -1 = flipped row order
        int start_y;            /// Beginning Y
        int end_y;              /// Ending Y
        int offbits;            /// Offset to image data
        uchar bit;              /// Bits in image
        uchar byte;             /// Bytes in image
        uchar*ptr;              /// Pointer into pixels
        uchar colormap[256][3]; /// Colormap
        uchar havemask = 0;     /// Single bit mask follows image data
        int use_5_6_5 = 0;      /// Use 5:6:5 for R:G:B channels in 16 bit images
        uint32_t w = 0;
        uint32_t h = 0;
        
        // Read offset to image data
        memcpy( &offbits, &buffer[buffque], 4 );
        buffque += 4;
        
        // Then the bitmap information...
        memcpy( &info_size, &buffer[buffque], 4 );
        buffque += 4;
        
        if (info_size < 40) 
        {
            // Old Windows/OS2 BMP header...
            w = flimgtk_memread_word( &buffer[buffque], &buffque );
            h = flimgtk_memread_word( &buffer[buffque], &buffque );
            
            // skip a WORD
            buffque += 2;

            depth = flimgtk_memread_word( &buffer[buffque], &buffque );
            
            compression = FLIMGTK_BI_RGB;
            colors_used = 0;

            repcount = info_size - 12;
        } 
        else 
        {
            // New BMP header...
            w = flimgtk_memread_dword( &buffer[buffque], &buffque );
            
            // If the height is negative, the row order is flipped
            temp = flimgtk_memread_int( &buffer[buffque], &buffque );
            
            if (temp < 0)
            {
                row_order = 1;
            }
            
            h = abs(temp);
            
            // Skip a WORD
            buffque += 2;
            
            depth = flimgtk_memread_word( &buffer[buffque], &buffque );
            compression = flimgtk_memread_dword( &buffer[buffque], &buffque );
            dataSize = flimgtk_memread_dword( &buffer[buffque], &buffque );
            
            // Skip a couple of DWORD
            buffque += 8;

            colors_used = flimgtk_memread_dword( &buffer[buffque], &buffque );
            
            // Skip DWORD
            buffque += 4;

            repcount = info_size - 40;

            if (!compression && depth>=8 && w>32/depth) 
            {
                int Bpp = depth/8;
                int maskSize = (((w*Bpp+3)&~3)*h) + (((((w+7)/8)+3)&~3)*h);
                if (maskSize==2*dataSize) 
                {
                    havemask = 1;
                    h /= 2;
                    bDepth = 4;
                }
            }
        }
        
        // Buffer skip by repcount.
        if ( repcount > 0 )
        {
            buffque += repcount;
            repcount = 0;
        }

        // check w,h,d
        if ( ( w == 0 ) || ( h == 0 ) || ( depth == 0 ) )
        {
            return NULL;
        }
        
        // Color map ?
        if (colors_used == 0 && depth <= 8)
        {
            colors_used = 1 << depth;
        }
        
        for (repcount = 0; repcount < colors_used; repcount ++) 
        {
            // Read BGR color palette...
            memcpy( &colormap[repcount], &buffer[buffque], 3 );
            buffque += 3;
            
            // Skip pad byte for new BMP files...
            if (info_size > 12)
            {
                buffque++;
            }
        }
        
        // Read first dword of colormap. It tells us if 5:5:5 or 5:6:5 for 16 bit
        if (depth == 16)
        {
            uint32_t tmpdw = flimgtk_memread_dword( &buffer[buffque], &buffque );
            if ( tmpdw == 0x0000f800 )
            {
                use_5_6_5 = 1;
            }
        }

        // Set byte depth for RGBA images
        if ( depth == 32 )
        {
            bDepth = 4;
        }

        // Setup image and buffers...
        if ( offbits > 0 ) 
        {
            buffque = offbits;
        }

        //---------------
        
        uchar* array = new uchar[ w * h * bDepth ];
        
        if ( array == NULL )
        {
            return NULL;
        }
        
        // Read the image data...
        color    = 0;
        repcount = 0;
        align    = 0;
        byte     = 0;
        temp     = 0;

        if (row_order < 0) 
        {
            start_y = h - 1;
            end_y   = -1;
        } 
        else 
        {
            start_y = 0;
            end_y   = h;
        }

        for (y = start_y; y != end_y; y += row_order) 
        {
            ptr = (uchar*)array + y * w * bDepth;

            switch( depth )
            {
                case 1 : /// Bitmap
                    for ( x = w, bit = 128; x>0; x-- ) 
                    {
                        if ( bit == 128 ) 
                        {
                            byte = buffer[buffque++];
                        }

                        if (byte & bit) 
                        {
                            *ptr++ = colormap[1][2];
                            *ptr++ = colormap[1][1];
                            *ptr++ = colormap[1][0];
                        } 
                        else 
                        {
                            *ptr++ = colormap[0][2];
                            *ptr++ = colormap[0][1];
                            *ptr++ = colormap[0][0];
                        }

                        if ( bit > 1 )
                        {
                            bit >>= 1;
                        }
                        else
                        {
                            bit = 128;
                        }
                    }

                    // Read remaining bytes to align to 32 bits...
                    for ( temp = ( w + 7 ) / 8; temp & 3; temp++ )
                    {
                        buffque++;
                    }
                    
                    break;

                case 4 : /// 16-color
                    for ( x = w, bit = 0xf0; x>0; x-- ) 
                    {
                        // Get a new repcount as needed...
                        if (repcount == 0) 
                        {
                            if (compression != FLIMGTK_BI_RLE4) 
                            {
                                repcount = 2;
                                color = -1;
                            } 
                            else 
                            {
                                while (align > 0) 
                                {
                                    align --;
                                    buffque++;
                                }

                                if ( (repcount = buffer[buffque++]) == 0) 
                                {
                                    if ( (repcount = buffer[buffque++]) == 0) 
                                    {
                                        // End of line...
                                            x ++;
                                        continue;
                                    } 
                                    else if (repcount == 1) 
                                    {
                                        // End of image...
                                        break;
                                    } 
                                    else if (repcount == 2) 
                                    {
                                        // Delta...
                                        repcount = buffer[buffque] *
                                                   buffer[buffque+1] *
                                                   w;
                                        buffque += 2;
                                        color = 0;
                                    } 
                                    else 
                                    {
                                        // Absolute...
                                        color = -1;
                                        align = ((4 - (repcount & 3)) / 2) & 1;
                                    }
                                } 
                                else 
                                {
                                    color = buffer[buffque++];
                                }
                            }
                        }

                        // Get a new color as needed...
                        repcount --;

                        // Extract the next pixel...
                        if (bit == 0xf0) 
                        {
                            // Get the next color byte as needed...
                            if (color < 0) 
                            {
                                temp = (uchar)buffer[buffque++];
                            }
                            else
                            {
                                temp  = (uchar)color;
                            }

                            // Copy the color value...
                            *ptr++ = colormap[(temp >> 4) & 15][2];
                            *ptr++ = colormap[(temp >> 4) & 15][1];
                            *ptr++ = colormap[(temp >> 4) & 15][0];

                            bit  = 0x0f;
                        } 
                        else 
                        {
                            bit  = 0xf0;

                            // Copy the color value...
                            *ptr++ = colormap[temp & 15][2];
                            *ptr++ = colormap[temp & 15][1];
                            *ptr++ = colormap[temp & 15][0];
                        }

                    }

                    if (!compression) 
                    {
                        // Read remaining bytes to align to 32 bits...
                        for (temp = (w + 1) / 2; temp & 3; temp ++) 
                        {
                            buffque++;
                        }
                    }
                    break;

                case 8 : /// 256-color
                    for ( x=w; x>0; x-- ) 
                    {
                        // Get a new repcount as needed...
                        if ( compression != FLIMGTK_BI_RLE8 ) 
                        {
                            repcount = 1;
                            color = -1;
                        }

                        if (repcount == 0) 
                        {
                            while (align > 0) 
                            {
                                align --;
                                buffque++;
                            }

                            if ((repcount = buffer[buffque++]) == 0) 
                            {
                                if ((repcount = buffer[buffque++]) == 0) 
                                {
                                    // End of line...
                                    x ++;
                                    continue;
                                } else if (repcount == 1) 
                                {
                                    // End of image...
                                    break;
                                } else if (repcount == 2) 
                                {
                                    // Delta...
                                    repcount = buffer[buffque] * 
                                               buffer[buffque+1] * 
                                               w;
                                    buffque += 2;
                                    color = 0;
                                } 
                                else 
                                {
                                    // Absolute...
                                    color = -1;
                                    align = (2 - (repcount & 1)) & 1;
                                }
                            } 
                            else 
                            {
                                color = buffer[buffque++];
                            }
                        }

                        // Get a new color as needed...
                        if (color < 0) 
                        {
                            temp = (uchar)buffer[buffque++];
                        }
                        else 
                        {
                            temp = (uchar)color;
                        }

                        // Check Error range ...
                        if ( ( temp < 0 ) || ( temp >= (int)colors_used ) )
                        {
                            temp = (int)colors_used - 1;
                        }

                        repcount --;

                        // Copy the color value...
                        *ptr++ = colormap[temp][2];
                        *ptr++ = colormap[temp][1];
                        *ptr++ = colormap[temp][0];
                        
                        if (havemask) 
                        {
                            ptr++;
                        }
                    }

                    if (!compression) 
                    {
                        // Read remaining bytes to align to 32 bits...
                        for ( temp=w; temp & 3; temp++ ) 
                        {
                            buffque++;
                        }
                    }
                    break;

                case 16 : /// 16-bit 5:5:5 or 5:6:5 RGB
                    for ( x=w; x>0; x--, ptr+=bDepth ) 
                    {
                        uchar b = buffer[buffque++];
                        uchar a = buffer[buffque++];
                        
                        if (use_5_6_5) 
                        {
                            ptr[2] = (uchar)(( b << 3 ) & 0xf8);
                            ptr[1] = (uchar)(((a << 5) & 0xe0) | ((b >> 3) & 0x1c));
                            ptr[0] = (uchar)(a & 0xf8);
                        } else {
                            ptr[2] = (uchar)((b << 3) & 0xf8);
                            ptr[1] = (uchar)(((a << 6) & 0xc0) | ((b >> 2) & 0x38));
                            ptr[0] = (uchar)((a<<1) & 0xf8);
                        }
                    }

                    // Read remaining bytes to align to 32 bits...
                    for ( temp=w * 2; temp & 3; temp ++ ) 
                    {
                        buffque++;
                    }
                    break;

                case 24 : /// 24-bit RGB
                    for ( x=w; x>0; x--, ptr += bDepth ) 
                    {
                        ptr[2] = (uchar)buffer[buffque++];
                        ptr[1] = (uchar)buffer[buffque++];
                        ptr[0] = (uchar)buffer[buffque++];
                    }

                    // Read remaining bytes to align to 32 bits...
                    for (temp = w * 3; temp & 3; temp ++) 
                    {
                        buffque++;
                    }
                    break;

                case 32 : /// 32-bit RGBA
                    for ( x=w; x>0; x--, ptr += bDepth ) 
                    {
                        ptr[2] = (uchar)buffer[buffque++];
                        ptr[1] = (uchar)buffer[buffque++];
                        ptr[0] = (uchar)buffer[buffque++];
                        ptr[3] = (uchar)buffer[buffque++];
                    }
                    break;
            } /// of switch(depth) ...
        } /// of for() ...
  
        if (havemask) 
        {
            for (y = h; y-- != 0; ) 
            {
                ptr = (uchar *)array + y * w * bDepth + 3;
                
                for( x = (w + 1), bit = 128; x-- != 0; ptr+=bDepth ) 
                {
                    if (bit == 128) 
                    {
                        byte = (uchar)buffer[buffque++];
                    }
                    
                    if (byte & bit)
                    {
                        *ptr = 0;
                    }
                    else
                    {
                        *ptr = 255;
                    }
                    
                    if (bit > 1)
                    {
                        bit >>= 1;
                    }
                    else
                    {
                        bit = 128;
                    }
                }
                // Read remaining bytes to align to 32 bits...
                for (temp = (w + 7) / 8; temp & 3; temp ++)
                {
                    buffque++;
                }
            }
        }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( array, w, h, bDepth );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else        
        return new Fl_RGB_Image( array, w, h, bDepth );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }
        
    return NULL;
}

Fl_RGB_Image* fl_imgtk::fliphorizontal( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    OMPSIZE_T w = (size_t)img->w();
    OMPSIZE_T h = (size_t)img->h();
    OMPSIZE_T d = (size_t)img->d();

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        const uchar* ptr = (const uchar*)img->data()[0];
        uchar* buff = new uchar[ w * h * d ];

        if ( buff == NULL )
            return NULL;

        memcpy( buff, ptr, w * h * d );

        OMPSIZE_T hcenter = h/2;
        OMPSIZE_T cnth = 0;
        OMPSIZE_T cntw = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<w; cntw++ )
        {
            for( cnth=0; cnth<hcenter; cnth++ )
            {
                size_t pos1 = ( w * ( h - 1 - cnth ) + cntw ) * d;
                size_t pos2 = ( w * cnth + cntw ) * d;

                fl_imgtk_swap_mem( &buff[ pos1 ], &buff[ pos2 ], d );
            }
        }        
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buff, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

bool fl_imgtk::fliphorizontal_ex( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return false;

    uchar* ptr = (uchar*)img->data()[0];
    OMPSIZE_T w = img->w();
    OMPSIZE_T h = img->h();
    OMPSIZE_T d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        OMPSIZE_T hcenter = h/2;
        OMPSIZE_T cnth = 0;
        OMPSIZE_T cntw = 0;

        #pragma omp parallel for private(cntw)
        for( cntw=0; cntw<w; cntw++ )
        {
            for( cnth=0; cnth<hcenter; cnth++ )
            {
                unsigned pos1 = ( w * ( h - 1 - cnth ) + cntw ) * d;
                unsigned pos2 = ( w * cnth + cntw ) * d;

                fl_imgtk_swap_mem( &ptr[ pos1 ], &ptr[ pos2 ], d );
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

    OMPSIZE_T w = img->w();
    OMPSIZE_T h = img->h();
    OMPSIZE_T d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        const uchar* ptr = (const uchar*)img->data()[0];
        uchar* buff = new uchar[ w * h * d ];

        if ( buff == NULL )
            return NULL;

        memcpy( buff, ptr, w * h * d );

        OMPSIZE_T wcenter = w/2;
        OMPSIZE_T cntw = 0;
        OMPSIZE_T cnth = 0;

        #pragma omp parallel for private(cntw)
        for( cnth=0; cnth<h; cnth++ )
        {
            for( cntw=0; cntw<wcenter; cntw++ )
            {
                unsigned pos1 = ( w * cnth + cntw ) * d;
                unsigned pos2 = ( w * cnth + ( w - cntw - 1 ) ) * d;

                fl_imgtk_swap_mem( &buff[ pos1 ], &buff[ pos2 ], d );
            }
        }

#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buff, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

bool fl_imgtk::flipvertical_ex( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return false;

    OMPSIZE_T w = img->w();
    OMPSIZE_T h = img->h();
    OMPSIZE_T d = img->d();

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        uchar* ptr = (uchar*)img->data()[0];
        OMPSIZE_T wcenter = w/2;
        OMPSIZE_T cntw = 0;
        OMPSIZE_T cnth = 0;

        #pragma omp parallel for private(cntw)
        for( cnth=0; cnth<h; cnth++ )
        {
            for( cntw=0; cntw<wcenter; cntw++ )
            {
                unsigned pos1 = ( w * cnth + cntw ) * d;
                unsigned pos2 = ( w * cnth + ( w - cntw - 1 ) ) * d;

                fl_imgtk_swap_mem( &ptr[ pos1 ], &ptr[ pos2 ], d );
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

    OMPSIZE_T w = img->w();
    OMPSIZE_T h = img->h();
    OMPSIZE_T d = img->d();

    OMPSIZE_T src_w = w;
    OMPSIZE_T src_h = h;

    if ( ( src_w > 0 ) && ( src_h > 0 ) )
    {
        const uchar* ptr = (const uchar*)img->data()[0];
        OMPSIZE_T new_w = src_h;
        OMPSIZE_T new_h = src_w;

        uchar* buff = new uchar[ new_w * new_h * d ];

        if ( buff == NULL )
            return NULL;

        OMPSIZE_T cntw = 0;
        OMPSIZE_T cnth = 0;

        for( cntw=new_w; cntw-- != 0; )
        {
            #pragma omp parallel for
            for( cnth=0; cnth<new_h; cnth++ )
            {
                unsigned pos1 = ( new_w * cnth + cntw ) * d;
                unsigned pos2 = ( src_w * ( new_w - cntw - 1 ) + cnth ) * d;

                memcpy( &buff[ pos1 ], &ptr[ pos2 ], d );
            }
        }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, new_w, new_h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buff, new_w, new_h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::rotate180( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    OMPSIZE_T w = img->w();
    OMPSIZE_T h = img->h();
    OMPSIZE_T d = img->d();

    OMPSIZE_T cur_w = w;
    OMPSIZE_T cur_h = h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        uchar* ptr = (uchar*)img->data()[0];
        uchar* buff = new uchar[ w * h * d ];

        if ( buff == NULL )
            return NULL;

        memcpy( buff, ptr, w * h * d );

        OMPSIZE_T imgmax = w*h;
        OMPSIZE_T cntmax = imgmax / 2;
        
        #pragma omp parallel for
        for( OMPSIZE_T cnt=0; cnt<cntmax; cnt++ )
        {
            fl_imgtk_swap_mem( &buff[ cnt * d ],
                               &buff[ (imgmax - cnt - 1) * d ],
                               d );
        }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buff, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::rotate270( Fl_RGB_Image* img )
{
    if ( img == NULL )
        return NULL;

    OMPSIZE_T w = img->w();
    OMPSIZE_T h = img->h();
    OMPSIZE_T d = img->d();

    OMPSIZE_T src_w = w;
    OMPSIZE_T src_h = h;

    if ( ( src_w > 0 ) && ( src_h > 0 ) )
    {
        const uchar* ptr = (const uchar*)img->data()[0];
        OMPSIZE_T new_w = src_h;
        OMPSIZE_T new_h = src_w;

        uchar* buff = new uchar[ new_w * new_h * d ];

        if ( buff == NULL )
            return NULL;

        OMPSIZE_T cntw = 0;
        OMPSIZE_T cnth = 0;

        #pragma omp parallel for private( cnth )
        for( cntw=0; cntw<new_w; cntw++ )
        {
            for( cnth=new_h; cnth-- != 0; )
            {
                OMPSIZE_T pos1 = ( new_w * cnth + cntw ) * d;
                OMPSIZE_T pos2 = ( src_w * cntw + new_h - cnth - 1 ) * d;

                memcpy( &buff[ pos1 ], &ptr[ pos2 ], d );
            }
        }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, new_w, new_h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buff, new_w, new_h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

inline float fl_imgtk_min4(float a, float b, float c, float d)
{
   float mn = a;
   if(mn > b) mn = b;
   if(mn > c) mn = c;
   if(mn > d) mn = d;
   return mn;
}

inline float fl_imgtk_max4(float a, float b, float c, float d)
{
   float mx = a;
   if(mx < b) mx = b;
   if(mx < c) mx = c;
   if(mx < d) mx = d;
   return mx;
}

// rotatefree() Code inspired from
// http://www.codeguru.com/cpp/g-m/gdi/article.php/c3693/ ...
// Rotate-a-Bitmap-at-Any-Angle-Without-GetPixelSetPixel.htm
// Maybe, it will slowly works if not enable openMP.
Fl_RGB_Image* fl_imgtk::rotatefree( Fl_RGB_Image* img, float deg )
{
    // I'm not sure this will be work for depth 1 or 2 ...
    // not do not test depth for a monment until it fully tested.
    if ( img == NULL )
        return NULL;

    long img_w = img->w();
    long img_h = img->h();
    long img_d = img->d();

    float CtX = ( (float)img_w ) / 2.0f;
    float CtY = ( (float)img_h ) / 2.0f;

    float fdeg = fl_imgtk_degree2f( deg );

    float cA = (float)cos( fdeg );
    float sA = (float)sin( fdeg );

    float x1 = CtX + (-CtX) * cA - (-CtY) * sA;
    float x2 = CtX + ((float)img_w - CtX) * cA - (-CtY) * sA;
    float x3 = CtX + ((float)img_w - CtX) * cA - ((float)img_h - CtY) * sA;
    float x4 = CtX + (-CtX) * cA - (img_h - CtY) * sA;

    float y1 = CtY + (-CtY) * cA + (-CtX) * sA;
    float y2 = CtY + ((float)img_h - CtY) * cA + (-CtX) * sA;
    float y3 = CtY + ((float)img_h - CtY) * cA + ((float)img_w - CtX) * sA;
    float y4 = CtY + (-CtY) * cA + ((float)img_w - CtX) * sA;

    long OfX = (long)floor(fl_imgtk_min4(x1, x2, x3, x4));
    long OfY = (long)floor(fl_imgtk_min4(y1, y2, y3, y4));

    long dstW = (long)ceil(fl_imgtk_max4(x1, x2, x3, x4)) - OfX;
    long dstH = (long)ceil(fl_imgtk_max4(y1, y2, y3, y4)) - OfY;

    // Now new image !
    uchar* obuff = new uchar[ dstW * dstH * img_d ];

    if ( obuff == NULL )
        return NULL;

    memset( obuff, 0x00, dstW * dstH * img_d );

    uchar* psrc = (uchar*)img->data()[0];

    long stepY = 0;
    long stepX = 0;

    // pointer to destination.
    uchar* dst = obuff;
    
    #pragma omp parallel for private( stepX )
    for ( stepY=0; stepY<dstH; stepY++ )
    {
        for ( stepX=0; stepX<dstW; stepX++ )
        {
#if defined(USING_INTERPOLATED_ROTATE_FREE)
            float CtX2 = CtX - OfX;
            float CtY2 = CtY - OfY;

            float orgX = ( cA*(stepX-CtX2) + sA*(stepY-CtY2)) + CtX;
            float orgY = (-sA*(stepX-CtX2) + cA*(stepY-CtY2)) + CtY;

            long iorgX  = (long) orgX;
            long iorgY  = (long) orgY;

            float diffX = (orgX - (float)iorgX);
            float diffY = (orgY - (float)iorgY);

            if ( ( (orgX >= 0) && (orgY >= 0) ) && \
                 ( (orgX < img_w-1) && (orgY < img_h-1) ) )
            {
                uchar* pd = &obuff[ ( stepY * dstW + stepX ) * img_d ];
                uchar* ps = &psrc[ ( iorgY * img_w + iorgX ) * img_d ];

                // Doing interpolated pixel calculation .
                for( unsigned cntd=0; cntd<img_d; cntd++ )
                {
                    float pv[4] = {0.0f};

                    // steps in order :
                    //   0 : current position.
                    //   1 : right side.
                    //   2 : below side.
                    //   3 : right below side.
                    pv[0] = (float)ps[ cntd ];
                    pv[1] = (float)ps[ cntd + img_d ];
                    pv[2] = (float)ps[ cntd + ( img_w * img_d ) ];
                    pv[3] = (float)ps[ cntd + ( ( img_w + 1 ) * img_d ) ];

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

            if ( ( (iorgX >= 0) && (iorgY >= 0) ) \
                 && ( (iorgX < img_w) && (iorgY < img_h) ) )
            {
                memcpy( &dst[ stepX * img_d ],
                        &psrc[ ( iorgX + iorgY * img_w ) * img_d ],
                        img_d );
            }
#endif
        }

        dst += dstW * img_d;
    }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    Fl_RGB_Image* newimg = new Fl_RGB_Image( obuff, dstW, dstH, img_d );
    if ( newimg != NULL )
    {
        newimg->alloc_array = 1;
        return newimg;
    }
#else
    return new Fl_RGB_Image( obuff, dstW, dstH, img->d() );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)

    // prevent compiler warning.
    return NULL;
}

Fl_RGB_Image* fl_imgtk_curve( Fl_RGB_Image* img, const uchar* LUT )
{
    if ( img == NULL )
        return NULL;

    const uchar* ptr = (const uchar*)img->data()[0];
    OMPSIZE_T w  = img->w();
    OMPSIZE_T h  = img->h();
    OMPSIZE_T d  = img->d();
    OMPSIZE_T td = d;
    OMPSIZE_T imgsz = w*h;

    if ( imgsz == 0 )
        return NULL;

    uchar* buff = new uchar[ imgsz * d ];

    if ( buff == NULL )
        return NULL;

    memcpy( buff, ptr, imgsz * d );

    if ( ( d == 2 ) || ( d == 4 ) ) td--;
    
    #pragma omp parallel for
    for( OMPSIZE_T cnt=0; cnt<imgsz; cnt++ )
    {
        for( OMPSIZE_T cntd=0; cntd<td; cntd++ )
        {
            buff[ cnt * d + cntd ] = LUT[ buff[ cnt * d + cntd ] ];
        }
    }
    
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    Fl_RGB_Image* newimg = new Fl_RGB_Image( buff, w, h, d );
    if ( newimg != NULL )
    {
        newimg->alloc_array = 1;
        return newimg;
    }
#else
    return new Fl_RGB_Image( buff, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)

    // prevent compiler warning.
    return NULL;
}

bool fl_imgtk_curve_ex( Fl_RGB_Image* img, const uchar* LUT )
{
    if ( img == NULL )
        return false;

    uchar* ptr   = (uchar*)img->data()[0];
    OMPSIZE_T w  = img->w();
    OMPSIZE_T h  = img->h();
    OMPSIZE_T d  = img->d();
    OMPSIZE_T td = d;
    OMPSIZE_T imgsz = w*h;

    if ( imgsz == 0 )
        return false;

    // sense alpha channel ( alpha channels will be skipped )
    if (( d == 2 ) || ( d == 4 )) td--;
    
    #pragma omp parallel for
    for( OMPSIZE_T cnt=0; cnt<imgsz; cnt++ )
    {
        for( OMPSIZE_T cntd=1; cntd<=td; cntd++ )
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

        lut[ cnt ] = (uchar)floor( col + 0.5f );
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

bool fl_imgtk::brightness_ex( Fl_RGB_Image* img, double perc )
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
        OMPSIZE_T img_w = img->w();
        OMPSIZE_T img_h = img->h();
        OMPSIZE_T img_d = img->d();
        OMPSIZE_T td    = img_d;

        if ( ( img_w == 0 ) || ( img_h == 0 ) )
            return NULL;
            
        // skip alpha channel.
        if ( ( img_d == 2 ) || ( img_d == 4 ) ) 
            td--;

        OMPSIZE_T buffsz = img_w * img_h;
        uchar* newbuff = new uchar[ buffsz * img_d ];

        if ( newbuff != NULL )
        {
            uchar* refbuff = (uchar*)img->data()[0];

            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<buffsz; cnt++ )
            {
                uchar* psrc = &refbuff[ cnt * img_d ];
                uchar* pdst = &newbuff[ cnt * img_d ];

                // alpha channel will skipped to invert.
                for( OMPSIZE_T rpt=0; rpt<td; rpt++ )
                {
                    pdst[ rpt ] = 0xFF - psrc[ rpt ];
                }
                
                if ( img_d != td )
                {
                    pdst[ td ] = psrc[ td ];
                }
            }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            Fl_RGB_Image* newimg = new Fl_RGB_Image( newbuff, img_w, img_h, img_d );
            if ( newimg != NULL )
            {
                newimg->alloc_array = 1;
                return newimg;
            }
#else
            return new Fl_RGB_Image( newbuff, img_w, img_h, img_d );
#endif /// of  #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        }
    }

    return NULL;
}

bool fl_imgtk::invert_ex( Fl_RGB_Image* img )
{
    if ( img != NULL )
    {
        OMPSIZE_T img_w = img->w();
        OMPSIZE_T img_h = img->h();
        OMPSIZE_T img_d = img->d();
        OMPSIZE_T td    = img_d;

        if ( ( img_w == 0 ) || ( img_h == 0 ) )
            return false;
            
        // skip alpha channel
        if ( ( img_d == 2 ) || ( img_d == 4 ) ) 
            td--;

        OMPSIZE_T buffsz = img_w * img_h;

        uchar* refbuff = (uchar*)img->data()[0];

        #pragma omp parallel for
        for( OMPSIZE_T cnt=0; cnt<buffsz; cnt++ )
        {
            uchar* psrc = &refbuff[ cnt * img_d ];

            // alpha channel will skipped to invert.
            for( OMPSIZE_T rpt=0; rpt<td; rpt++ )
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
        OMPSIZE_T img_w = img->w();
        OMPSIZE_T img_h = img->h();
        OMPSIZE_T img_d = img->d();
        OMPSIZE_T td    = img_d;

        if ( ( img_w == 0 ) || ( img_h == 0 ) )
            return NULL;

        if ( ( kfc->w == 0 ) || ( kfc->h == 0 ) || ( kfc->msz == 0 ) || ( kfc->m == NULL ) )
            return NULL;

        uchar* pixels  = (uchar*)img->data()[0];
        uchar* newbuff = new uchar[ img_w * img_h * img_d ];

        if ( newbuff == NULL )
            return NULL;

        if ( ( img_d == 2 ) || ( img_d == 4 ) )
            td--;

        #pragma omp parallel for
        for( OMPSIZE_T cntx=0; cntx<img_w; cntx++ )
        {
            for( OMPSIZE_T cnty=0; cnty<img_h; cnty++ )
            {
                double adj[4] = {0.0};

                // -- applying matrix ---
                for( OMPSIZE_T fcntx=0; fcntx<kfc->w; fcntx++ )
                {
                    for( OMPSIZE_T fcnty=0; fcnty<kfc->h; fcnty++ )
                    {
                        OMPSIZE_T posX = ( cntx - kfc->w / 2 + fcntx + img_w )
                                         % img_w;
                        OMPSIZE_T posY = ( cnty - kfc->h / 2 + fcnty + img_h )
                                         % img_h;

                        OMPSIZE_T posM = posY * img_w + posX;

                        if ( posM < ( img_w * img_h ) )
                        {
                            for( OMPSIZE_T cntd=0; cntd<td; cntd ++ )
                            {
                                adj[ cntd ] += (double)pixels[ posM * img_d  + cntd ] *
                                               (double)kfc->m[ fcnty * kfc->w + fcntx ];
                            }
                        }
                    }
                }

                for( OMPSIZE_T cntd=0; cntd<td; cntd++ )
                {
                    uchar rpixel = MIN( MAX( kfc->f * adj[ cntd ] + kfc->b, 0 ), 255 );
                    newbuff[ ( cnty * img_w + cntx ) * img_d + cntd ] = rpixel;
                }
                
                if ( img_d != td ) /// alpha channel
                {
                    size_t pq = ( cnty * img_w + cntx ) * img_d + td;
                    newbuff[ pq ] = pixels[ pq ]; 
                }
            }
        }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( newbuff, img_w, img_h, img_d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( newbuff, img_w, img_h, img_d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

bool fl_imgtk::filtered_ex( Fl_RGB_Image* img, kfconfig* kfc )
{
    if ( ( img != NULL ) && ( kfc != NULL ) )
    {
        OMPSIZE_T img_w = img->w();
        OMPSIZE_T img_h = img->h();
        OMPSIZE_T img_d = img->d();
        OMPSIZE_T td    = img_d;

        if ( ( img_w == 0 ) || ( img_h == 0 ) )
            return false;

        if ( ( kfc->w == 0 ) || ( kfc->h == 0 ) || ( kfc->msz == 0 ) || ( kfc->m == NULL ) )
            return false;
            
        if ( ( img_d == 2 ) || ( img_d == 4 ) )
            td--;

        uchar* pixels  = (uchar*)img->data()[0];

        #pragma omp parallel for
        for( OMPSIZE_T cntx=0; cntx<img_w; cntx++ )
        {
            for( OMPSIZE_T cnty=0; cnty<img_h; cnty++ )
            {
                double adj[4] = {0.0};

                // -- applying matrix ---
                for( OMPSIZE_T fcntx=0; fcntx<kfc->w; fcntx++ )
                {
                    for( OMPSIZE_T fcnty=0; fcnty<kfc->h; fcnty++ )
                    {
                        OMPSIZE_T posX = ( cntx - kfc->w / 2 + fcntx + img_w )
                                         % img_w;
                        OMPSIZE_T posY = ( cnty - kfc->h / 2 + fcnty + img_h )
                                         % img_h;

                        OMPSIZE_T posM = posY * img_w + posX;

                        if ( posM < ( img_w * img_h ) )
                        {
                            for( OMPSIZE_T cntd=0; cntd<td; cntd++ )
                            {
                                adj[ cntd ] += (double)pixels[ posM * img_d  + cntd ] *
                                               (double)kfc->m[ fcnty * kfc->w + fcntx ];
                            }
                        }
                    }
                }

                for( OMPSIZE_T cntd=0; cntd<td; cntd++ )
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

bool fl_imgtk_gen_lowfreq( uchar** out, const uchar* src, unsigned w, unsigned h, unsigned d, unsigned sz )
{
    if ( src == NULL )
        return false;

    if ( ( w < sz*2 ) || ( h < sz*2 ) )
        return false;

    if ( sz < 2 )
        return false;

    OMPSIZE_T startpos = sz / 2;
    OMPSIZE_T endposw  = w - startpos;
    OMPSIZE_T endposh  = h - startpos;

    uchar* outbuff = new uchar[ w * h * d ];

    if ( outbuff == NULL )
        return false;

    OMPSIZE_T cnty;
    OMPSIZE_T cntx;
    OMPSIZE_T td = d;
    
    // skip alpha channel.
    if ( ( d == 2 ) || ( d == 4 ) ) td--;

    #pragma omp parallel for private( cntx )
    for( cnty=startpos; cnty<endposh; cnty++ )
    {
        for( cntx=startpos; cntx<endposw; cntx++ )
        {
            unsigned sum[4] = {0};

            for( OMPSIZE_T ry = -startpos; ry<(startpos+1); ry++ )
            {
                for( OMPSIZE_T rx = -startpos; rx<(startpos+1); rx++ )
                {
                    for( OMPSIZE_T rpt=0; rpt<td; rpt++ )
                    {
                        sum[rpt] += src[ ( ( cnty + ry ) * w + ( cntx + rx ) ) * d + rpt ];
                    }
                }
            }

            for( OMPSIZE_T rpt=0; rpt<td; rpt++ )
            {
                outbuff[ ( cnty * w + cntx ) * d + rpt ] \
                 = MIN( 255, sum[rpt] / ( sz * sz ) );
            }
            
            if ( td != d )
            {
                size_t bq = ( cnty * w + cntx ) * d + td;
                outbuff[ bq ] = src[ bq ];
            }
        }
    }

    *out = outbuff;

    return true;
}

Fl_RGB_Image* fl_imgtk::edgeenhance( Fl_RGB_Image* img, unsigned factor, unsigned margin )
{
    if ( img != NULL )
    {
        OMPSIZE_T imgsz = img->w() * img->h() * img->d();
        uchar* outbuff = new uchar[ imgsz ];

        if ( outbuff == NULL )
            return NULL;

        const uchar* rbuff  = (const uchar*)img->data()[0];
        uchar* lfimg5 = NULL;
        uchar* lfimg9 = NULL;

        if ( fl_imgtk_gen_lowfreq( &lfimg5, rbuff,
                                   img->w(), img->h(), img->d(), 5 ) == false )
            return NULL;

        if ( fl_imgtk_gen_lowfreq( &lfimg9, rbuff,
                                   img->w(), img->h(), img->d(), 9 ) == false )
        {
            delete[] lfimg5;
            return NULL;
        }

        float fedgev = (float)factor / 8.0f;

        if ( ( img->w() < ( margin * 2 ) ) || ( img->h() < ( margin * 2 ) ) )
        {
            margin = 0;
        }

        long cnty;
        long cntx;
        long cntw = img->w();
        long cnth = img->h();
        long mgnx = margin;
        long mgny = margin;
        long mgnw = img->w() - ( margin * 2 );
        long mgnh = img->h() - ( margin * 2 );
        long td   = (long)img->d();
        
        if ( ( img->d() == 2 ) || ( img->d() == 4 ) )
            td--;

        #pragma omp parallel for private( cntx )
        for( cnty=0; cnty<cnth; cnty++ )
        {
            for( cntx=0; cntx<cntw; cntx++ )
            {
                for( unsigned rpt=0; rpt<td; rpt++ )
                {
                    if ( ( cntx > mgnx ) && ( cntx < (mgnx+mgnw) ) && \
                         ( cnty > mgny ) && ( cnty < (mgny+mgnh) ) )
                    {
                        float dlv = (float)lfimg5[ ( cnty * img->w() + cntx ) * img->d() + rpt ]
                                    - (float)lfimg9[ ( cnty * img->w() + cntx ) * img->d() + rpt ];

                        unsigned pv = std::abs( (float)rbuff[ ( cnty * img->w() + cntx ) * img->d() + rpt ]
                                                + (float)dlv *  fedgev );

                        outbuff[ ( cnty * img->w() + cntx ) * img->d() + rpt ] = MIN( 255, pv );
                    }
                    else
                    {
                        outbuff[ ( cnty * img->w() + cntx ) * img->d() + rpt ] \
                        = rbuff[ ( cnty * img->w() + cntx ) * img->d() + rpt ];
                    }
                }
                
                if ( (long)img->d() != td ) /// alpha channel.
                {
                    long bq = ( cnty * img->w() + cntx ) * img->d() + td;
                    outbuff[ bq ] = rbuff[ bq ];
                }
            }
        }

        delete[] lfimg5;
        delete[] lfimg9;
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( outbuff, img->w(), img->h(), img->d() );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( outbuff, img->w(), img->h(), img->d() );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

bool fl_imgtk::edgeenhance_ex( Fl_RGB_Image* img, unsigned factor, unsigned margin )
{
    if ( img != NULL )
    {
        uchar* rbuff  = (uchar*)img->data()[0];
        uchar* lfimg5 = NULL;
        uchar* lfimg9 = NULL;

        if ( fl_imgtk_gen_lowfreq( &lfimg5, rbuff,
                                   img->w(), img->h(), img->d(), 5 ) == false )
            return false;

        if ( fl_imgtk_gen_lowfreq( &lfimg9, rbuff,
                                   img->w(), img->h(), img->d(), 9 ) == false )
        {
            delete[] lfimg5;
            return false;
        }

        float fedgev = (float)factor / 8.0f;

        long cnty;
        long cntx;
        long mgnx = margin;
        long mgny = margin;
        long mgnw = img->w() - ( margin * 2 );
        long mgnh = img->h() - ( margin * 2 );
        long td   = (long)img->d();

        #pragma omp parallel for private( cntx )
        for( cnty=mgny; cnty<mgnh; cnty++ )
        {
            for( cntx=mgnx; cntx<mgnw; cntx++ )
            {
                for( unsigned rpt=0; rpt<td; rpt++ )
                {
                    float dlv = (float)lfimg5[ ( cnty * img->w() + cntx ) * img->d() + rpt ]
                                - (float)lfimg9[ ( cnty * img->w() + cntx ) * img->d() + rpt ];

                    unsigned pv = std::abs( (float)rbuff[ ( cnty * img->w() + cntx ) * img->d() + rpt ]
                                            + (float)dlv *  fedgev );

                    rbuff[ ( cnty * img->w() + cntx ) * img->d() + rpt ] = MIN( 255, pv );
                }
            }
        }

        delete[] lfimg5;
        delete[] lfimg9;

        img->uncache();

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

Fl_RGB_Image* fl_imgtk::draw_currentwindow( void* w )
{
    Fl_Window* cwin = (Fl_Window*)w;
    
    if ( cwin == NULL )
    {
        cwin = Fl::first_window();
    }
    
    if ( cwin != NULL )
    {
        unsigned cwin_x = cwin->x();
        unsigned cwin_y = cwin->y();
        unsigned cwin_w = cwin->w();
        unsigned cwin_h = cwin->h();

        uchar* widgetbuff = new uchar[ cwin_w * cwin_h * 3 + 1 ];
        if ( widgetbuff != NULL )
        {
            // Read current window pixels to buffer for create a new Fl_RGB_Image.
            fl_read_image( widgetbuff, cwin_x, cwin_y,
                                       cwin_w, cwin_h );

            Fl_RGB_Image* newimg =  new Fl_RGB_Image( widgetbuff,
                                                      cwin_w,
                                                      cwin_h,
                                                      3 );
                                                      
            if ( newimg == NULL )
            {
                delete[] widgetbuff;
            }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            else
            {
                newimg->alloc_array = 1;
            }
#endif ///  of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            return newimg;
        }
    }
    
    return NULL;
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
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                        sdimg->alloc_array = 1;
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)                        
                        BSplineFilter* bsfilter = new BSplineFilter();
                        if ( bsfilter != NULL )
                        {
                            ResizeEngine* reup = new ResizeEngine( bsfilter );
                            if (  reup != NULL )
                            {
                                blurredimg = reup->scale( sdimg, w->w(), w->h() );
                                delete reup;
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                                if ( blurredimg != NULL )
                                {
                                    blurredimg->alloc_array = 1;
                                }
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
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
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                    sdimg->alloc_array = 1;
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                    BSplineFilter* bsfilter = new BSplineFilter();
                    if ( bsfilter != NULL )
                    {
                        ResizeEngine* reup = new ResizeEngine( bsfilter );
                        if (  reup != NULL )
                        {
                            newimg = reup->scale( sdimg, src->w(), src->h() );
                            delete reup;
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                            if ( newimg != NULL )
                            {
                                newimg->alloc_array = 1;
                            }
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
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
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                    sdimg->alloc_array = 1;
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                    BSplineFilter* bsfilter = new BSplineFilter();
                    if ( bsfilter != NULL )
                    {
                        ResizeEngine* reup = new ResizeEngine( bsfilter );
                        if (  reup != NULL )
                        {
                            Fl_RGB_Image* newimg = reup->scale( sdimg, src->w(), src->h() );

                            if ( newimg != NULL )
                            {
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
                                newimg->alloc_array = 1;
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
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
        OMPSIZE_T rsx = sx;
        OMPSIZE_T rsy = sy;
        OMPSIZE_T rw  = w;
        OMPSIZE_T rh  = h;
        OMPSIZE_T sd  = src->d();

        if ( ( rsx > src->w() ) || ( rsy > src->h() ) )
            return NULL;
        
        if ( ( ( rsx + w ) > src->w() ) || ( ( rsy + h ) > src->h() ) )
            return NULL;

        if ( src->w() < ( rw + rsx ) )
        {
            rw = src->w() - rsx;
        }

        if ( src->h() < ( rh + rsy ) )
        {
            rh = src->h() - rsy;
        }

        uchar* rbuff = (uchar*)src->data()[0];
        uchar* obuff = new uchar[ rw * rh * sd ];

        if ( obuff != NULL )
        {
            OMPSIZE_T srcw = rsx + w;
            OMPSIZE_T srch = rsy + h;
            OMPSIZE_T cnty;
            OMPSIZE_T cntx;
            OMPSIZE_T putx;
            OMPSIZE_T puty;

            // bug fixed thank you for, https://github.com/fire-eggs
            #pragma omp parallel for private( cntx )
            for( cnty=rsy; cnty<srch; cnty++ )
            {
                for( cntx=rsx; cntx<srcw; cntx ++ )
                {
                    putx = cntx - rsx;
                    puty = cnty - rsy;
                    uchar* rptr = &rbuff[ ( cnty * src->w() + cntx ) * sd ];
                    uchar* wptr = &obuff[ ( puty * rw + putx ) * sd ];

                    memcpy( wptr, rptr, sd );
                }
            }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
            Fl_RGB_Image* newimg = new Fl_RGB_Image( obuff, rw, rh , sd );
            if ( newimg != NULL )
            {
                newimg->alloc_array = 1;
                return newimg;
            }
#else
            return new Fl_RGB_Image( obuff, rw, rh , sd );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        }
    }

    return NULL;
}

bool fl_imgtk_convertrgb( uchar*& wbuff, unsigned w, unsigned h, unsigned d, uchar* rbuff, unsigned rd )
{
    if ( ( rbuff != NULL ) 
         && ( w > 0 ) && ( h > 0 ) && ( rd <= 2 ) && ( d >= 3 )
         && ( d != rd ) )
    {
        wbuff = new uchar[ w * h * d ];

        if ( wbuff == NULL )
            return false;

        if ( d == 3 )
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<(w*h); cnt++ )
            {
                wbuff[cnt*3+0] = \
                wbuff[cnt*3+1] = \
                wbuff[cnt*3+2] = rbuff[cnt];
            }
        }
        else
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<(w*h); cnt++ )
            {
                wbuff[cnt*4+0] = \
                wbuff[cnt*4+1] = \
                wbuff[cnt*4+2] = rbuff[cnt*2];
                wbuff[cnt*4+3] = rbuff[cnt*2+1];
            }            
        }
        
        return true;
    }
    return false;
}

// for effectible conversion ...
// all data will be converted to RGBA then converted to target depth.
bool fl_imgtk_putimgonbuffer( uchar* buff, unsigned bw, unsigned bh, unsigned bd,
                              Fl_RGB_Image* img, int px, int py, float alpha )
{
    if ( ( buff != NULL ) && ( img != NULL ) )
    {
        bool rbuffalloc    = false;
        bool wbuffalloc    = false;
        unsigned w_bd      = bd;
        uchar* wbuff       = buff;
        uchar* rbuff       = (uchar*)img->data()[0];
        Fl_RGB_Image* rimg = img;

        if ( w_bd < 3 ) // dest buffer convert to RGB
        {
            uchar* refb = buff;
        
            if ( w_bd == 1 ) /// no alpha
            {
                bool retb = \
                fl_imgtk_convertrgb( wbuff, bw, bh, 3, refb, w_bd );
                
                if ( retb == false )
                    return retb;

                w_bd = 3;
            }
            else /// has alpha
            {
                bool retb = \
                fl_imgtk_convertrgb( wbuff, bw, bh, 4, refb, w_bd );

                if ( retb == false )
                    return retb;

                w_bd = 4;
            }

            wbuffalloc = true;
        }

        if ( img->d() < 3 )
        {
            uchar* convbuff = NULL;
            int convdepth = 3;
            
            // convert to RGB or RGBA.
            if ( img->d() == 1 ) /// Y to RGB ( depth = 3 )
            {
                bool retb = \
                fl_imgtk_convertrgb( convbuff, 
                                     img->w(), img->h(), 3, 
                                     rbuff, 1 );
                                     
                if ( retb == false )
                {
                    if ( wbuffalloc == true )
                    {
                        delete[] wbuff;
                    }
                    return false;
                }
            }
            else /// Y+A to RGBA ( depth = 4 )
            {
                bool retb = \
                fl_imgtk_convertrgb( convbuff, 
                                     img->w(), img->h(), 4, 
                                     rbuff, 1 );
                                     
                if ( retb == false )
                {
                    if ( wbuffalloc == true )
                    {
                        delete[] wbuff;
                    }
                    return false;
                }                
                convdepth = 4;
            }
            
            // create converted image.
            rimg = new Fl_RGB_Image( convbuff, img->w(), img->h(), convdepth );
            if ( rimg == NULL )
            {
                delete[] convbuff;
                return false;
            }
            
            rimg->alloc_array = 1;
            rbuff = (uchar*)rimg->data()[0];
            rbuffalloc = true;
        }

        int minx = px;
        int miny = py;
        int imgx = 0;
        int imgy = 0;

        OMPSIZE_T maxw = MIN( rimg->w() + minx, bw );
        OMPSIZE_T maxh = MIN( rimg->h() + miny, bh );

        if ( px < 0 )
        {
            minx = 0;
            maxw = MIN( rimg->w() + px, bw );
            imgx = abs( px );
            px   = 0;
        }

        if ( py < 0 )
        {
            miny = 0;
            maxh = MIN ( rimg->h() + py, bh );
            imgy = abs( py );
            py   = 0;
        }

        OMPSIZE_T cntx = 0;
        OMPSIZE_T cnty = 0;
        OMPSIZE_T imgw = rimg->w();
        OMPSIZE_T imgd = rimg->d();
        OMPSIZE_T td   = imgd;
        
        if ( imgd == 4 ) /// let make skip alpha.
            td--;
            
        float a_opa = MIN( 1.f, alpha );
        
        #pragma omp parallel for private( cnty )
        for( cnty=miny; cnty<maxh; cnty++ )
        {
            for( cntx=minx; cntx<maxw; cntx++ )
            {
                unsigned ipx = imgx + (cntx - minx);
                unsigned ipy = imgy + (cnty - miny);

                uchar* rptr = &rbuff[ ( ipy * imgw + ipx ) * imgd ];
                uchar* wptr = &wbuff[ ( ( cnty * bw ) + cntx ) * w_bd ];

                // read-pixels
                float r_p[3] = {0.f};
                float r_a    = 1.0f;

                // write-pixels
                float w_p[3] = {0.f};
                float w_a    = 1.0f;
                
                for ( OMPSIZE_T dd=0; dd<3; dd++ )
                {
                    w_p[dd] = (float)wptr[dd] / 255.f;
                    r_p[dd] = (float)rptr[dd] / 255.f;
                }

                if ( w_bd == 4 )
                {
                    w_a = (float)wptr[3] / 255.f;
                }
                                                        
                if ( imgd == 4 )
                {
                    r_a = (float)rptr[3] / 255.f;
                }

                if ( a_opa < 1.f )
                {
                    r_a *= a_opa;
                }
                
                float  cTemp[3] = {0};
                
                for( OMPSIZE_T dd=0; dd<3; dd++ )
                {
                    cTemp[dd] = r_p[dd] * r_a + w_p[dd] * w_a * ( 1.0f - r_a );
                }
                
                r_a = w_a + ( 1.0f - w_a ) * r_a;

                if ( r_a == 0.0f )
                {
                    if ( a_opa > 0.0f )
                    {
                        for( OMPSIZE_T dd=0; dd<3; dd++ )
                        {    
                            cTemp[dd] = r_p[dd];
                        }
                    }
                    else
                    {
                        for( OMPSIZE_T dd=0; dd<3; dd++ )
                        {    
                            cTemp[dd] = w_p[dd];
                        }
                    }
                }
                else
                {
                    for( OMPSIZE_T dd=0; dd<3; dd++ )
                    {    
                        cTemp[dd] /= r_a;
                    }
                }
                
                if ( w_bd == 4 )
                {
                    for( OMPSIZE_T dd=0; dd<3; dd++ )
                    {
                        wptr[dd] = cTemp[dd] * 255.f;
                    }
                    wptr[3] = r_a * 255.f;
                }
                else
                {
                    for( OMPSIZE_T dd=0; dd<3; dd++ )
                    {
                        wptr[dd] = cTemp[dd] * ( 255.f * r_a );
                    }
                }
            }
        }
        
        if ( rbuffalloc == true )
        {
            delete rimg;
        }
        
        if ( wbuffalloc == true ) /// let convert to target.
        {
            // let convert to target depth.
            if ( bd == 1 ) /// depth 1 ( Y only )
            {
                #pragma omp parallel for
                for( OMPSIZE_T cnt=0; cnt<(bw*bh); cnt++ )
                {
                    ulong avrc = ( wbuff[cnt*3+0]
                                   + wbuff[cnt*3+1]
                                   + wbuff[cnt*3+2] ) / 3;
                    buff[cnt] = (uchar)(avrc & 0x000000FF);
                }
            }
            else /// must be depth 2 ( Y+A )
            {
                #pragma omp parallel for
                for( OMPSIZE_T cnt=0; cnt<(bw*bh); cnt++ )
                {
                    ulong avrc = ( wbuff[cnt*4+0]
                                   + wbuff[cnt*4+1]
                                   + wbuff[cnt*4+2] ) / 3;
                    buff[cnt*2] = (uchar)(avrc & 0x000000FF);
                    buff[cnt*2+1] = wbuff[cnt*4+3];
                }
            }
            
            delete[] wbuff;
        }
        
        return true;
    }
    
    return false;
}

bool fl_imgtk_subimgonbuffer( uchar* buff, unsigned bw, unsigned bh, unsigned bd,
                              Fl_RGB_Image* img, int px, int py, float alpha )
{
    if ( ( buff != NULL ) && ( img != NULL ) )
    {
        bool rbuffalloc    = false;
        bool wbuffalloc    = false;
        unsigned w_bd      = bd;
        uchar* wbuff       = buff;
        uchar* rbuff       = (uchar*)img->data()[0];
        Fl_RGB_Image* rimg = img;

        if ( w_bd < 3 ) // dest buffer convert to RGB
        {
            uchar* refb = buff;
            
            if ( w_bd == 1 ) /// no alpha
            {
                bool retb = \
                fl_imgtk_convertrgb( wbuff, bw, bh, 3, refb, w_bd );
                
                if ( retb == false )
                    return retb;

                w_bd = 3;
            }
            else /// has alpha
            {
                bool retb = \
                fl_imgtk_convertrgb( wbuff, bw, bh, 4, refb, w_bd );

                if ( retb == false )
                    return retb;

                w_bd = 4;
            }

            wbuffalloc = true;
        }

        if ( img->d() < 3 )
        {
            uchar* convbuff = NULL;
            int convdepth = 3;
            
            // convert to RGB or RGBA.
            if ( img->d() == 1 ) /// Y to RGB ( depth = 3 )
            {
                bool retb = \
                fl_imgtk_convertrgb( convbuff, 
                                     img->w(), img->h(), 3, 
                                     rbuff, 1 );
                                     
                if ( retb == false )
                {
                    if ( wbuffalloc == true )
                    {
                        delete[] wbuff;
                    }
                    return false;
                }
            }
            else /// Y+A to RGBA ( depth = 4 )
            {
                bool retb = \
                fl_imgtk_convertrgb( convbuff, 
                                     img->w(), img->h(), 4, 
                                     rbuff, 1 );
                                     
                if ( retb == false )
                {
                    if ( wbuffalloc == true )
                    {
                        delete[] wbuff;
                    }
                    return false;
                }                
                convdepth = 4;
            }
            
            // create converted image.
            rimg = new Fl_RGB_Image( convbuff, img->w(), img->h(), convdepth );
            if ( rimg == NULL )
            {
                delete[] convbuff;
                return false;
            }
            
            rimg->alloc_array = 1;
            rbuff = (uchar*)rimg->data()[0];
            rbuffalloc = true;
        }
        
        OMPSIZE_T minx = px;
        OMPSIZE_T miny = py;
        OMPSIZE_T imgx = 0;
        OMPSIZE_T imgy = 0;

        OMPSIZE_T maxw = MIN( rimg->w() + minx, bw );
        OMPSIZE_T maxh = MIN( rimg->h() + miny, bh );

        if ( px < 0 )
        {
            minx = 0;
            maxw = MIN( rimg->w() + px, bw );
            imgx = abs( px );
            px   = 0;
        }

        if ( py < 0 )
        {
            miny = 0;
            maxh = MIN ( rimg->h() + py, bh );
            imgy = abs( py );
            py   = 0;
        }

        OMPSIZE_T cntx = 0;
        OMPSIZE_T cnty = 0;
        OMPSIZE_T imgw = rimg->w();
        OMPSIZE_T imgd = rimg->d();

        #pragma omp parallel for private( cnty )
        for( cnty=miny; cnty<maxh; cnty++ )
        {
            for( cntx=minx; cntx<maxw; cntx++ )
            {
                OMPSIZE_T ipx = imgx + (cntx - minx);
                OMPSIZE_T ipy = imgy + (cnty - miny);

                uchar* rptr = &rbuff[ ( ipy * imgw + ipx ) * imgd ];
                uchar* wptr = &buff[ ( ( cnty * bw ) + cntx ) * bd ];

                uchar alp = 255;

                if ( imgd == 4 )
                {
                    alp = rptr[3];
                }

                float falp = ( (float)alp / 255.0f ) * alpha;

                for( OMPSIZE_T rpt=0; rpt<3; rpt++ )
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
            }
        } /// of for() --
        
        if ( rbuffalloc == true )
        {
            delete rimg;
        }
        
        if ( wbuffalloc == true ) /// let convert to target.
        {
            // let convert to target depth.
            if ( bd == 1 ) /// depth 1 ( Y only )
            {
                #pragma omp parallel for
                for( OMPSIZE_T cnt=0; cnt<(bw*bh); cnt++ )
                {
                    ulong avrc = ( wbuff[cnt*3+0]
                                   + wbuff[cnt*3+1]
                                   + wbuff[cnt*3+2] ) / 3;
                    buff[cnt] = (uchar)(avrc & 0x000000FF);
                }
            }
            else /// must be depth 2 ( Y+A )
            {
                #pragma omp parallel for
                for( OMPSIZE_T cnt=0; cnt<(bw*bh); cnt++ )
                {
                    ulong avrc = ( wbuff[cnt*4+0]
                                   + wbuff[cnt*4+1]
                                   + wbuff[cnt*4+2] ) / 3;
                    buff[cnt*2] = (uchar)(avrc & 0x000000FF);
                    buff[cnt*2+1] = wbuff[cnt*4+3];
                }
            }
            
            delete[] wbuff;
        }
        
        return true;
    }
    
    return false;
}

// result must be RGBA.
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
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
        }
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return newimg;
}

// result must be RGBA
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
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
        }
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
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

    OMPSIZE_T imgsz = src->w() * src->h();

    uchar* refbuff = (uchar*)src->data()[0];
    amap = new uchar[ imgsz ];

    if ( amap != NULL )
    {
        memset( amap, 0, imgsz );
        
        if ( src->d() == 4 )
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<imgsz; cnt++ )
            {
                amap[ cnt ] = (uchar)( (float)refbuff[cnt*4] * val );
            }
        }
        else
        {
            // Generate alphamap from RGB average.
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<imgsz; cnt++ )
            {
                float avrp = (float)refbuff[cnt*3+0] 
                             + (float)refbuff[cnt*3+1] 
                             + (float)refbuff[cnt*3+2];
                avrp /= 3.f;

                amap[ cnt ] = (uchar)( avrp * val );
            }
        }

        return imgsz;
    }

    return 0;
}

unsigned fl_imgtk::makealphamap( uchar* &amap, unsigned w, unsigned h, uchar val )
{
    unsigned imgsz = w * h;

    val = MAX( 1.0f, MIN( 0.0f, val ) );

    amap = new uchar[ imgsz ];

    if ( amap != NULL )
    {
        uchar fillv = (uchar)( 255.0f * val );
        memset( amap, fillv, imgsz );
        return imgsz;
    }

    return 0;
}

// it only availed for RGBA.
Fl_RGB_Image* fl_imgtk::applyalpha( Fl_RGB_Image* src, float val )
{
    Fl_RGB_Image* newimg = NULL;

    if ( src != NULL )
    {
        if ( src->d() < 3 )
            return NULL;

        OMPSIZE_T img_w  = src->w();
        OMPSIZE_T img_h  = src->h();
        OMPSIZE_T img_sz = img_w * img_h;

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
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* sbpos = &sbuff[ cnt * 3 ];
                uchar* obpos = &obuff[ cnt * 4 ];

                memcpy( obpos , sbpos, 3 );
                obpos[3] = 0xFF;
            }
        }

        #pragma omp parallel for
        for( OMPSIZE_T cnt=0; cnt<img_sz; cnt++ )
        {
            uchar* obp = &obuff[ cnt * 4 + 3 ];

            *obp = (uchar)( (float)*obp * val );
        }

        newimg = new Fl_RGB_Image( obuff, img_w, img_h, 4 );
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
        }
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return newimg;
}

// it only availed for RGBA.
Fl_RGB_Image* fl_imgtk::applyalpha( Fl_RGB_Image* src, uchar* alphamap, unsigned amsz )
{
    Fl_RGB_Image* newimg = NULL;

    if ( src != NULL )
    {
        if ( src->d() < 3 )
            return NULL;

        OMPSIZE_T img_w  = src->w();
        OMPSIZE_T img_h  = src->h();
        OMPSIZE_T img_sz = img_w * img_h;

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
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* sbpos = &sbuff[ cnt * 3 ];
                uchar* obpos = &obuff[ cnt * 4 ];

                memcpy( obpos , sbpos, 3 );
            }
        }

        // Apply alpha map.
        if ( ( alphamap != NULL ) && ( amsz == img_sz ) )
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<img_sz; cnt++ )
            {
                obuff[ cnt * 4 + 3 ] = alphamap[ cnt ];
            }
        }

        newimg = new Fl_RGB_Image( obuff, img_w, img_h, 4 );
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
        }
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return newimg;
}

// it only availed for RGBA.
bool fl_imgtk::applyalpha_ex( Fl_RGB_Image* src, float val )
{
    if ( src != NULL )
    {
        if ( src->d() < 3 )
            return false;

        OMPSIZE_T img_w  = src->w();
        OMPSIZE_T img_h  = src->h();
        OMPSIZE_T img_sz = img_w * img_h;

        uchar* ptr = (uchar*)src->data()[0];

        val = MIN( 1.0f, MAX( 0.0f, val ) );

        if ( src->d() == 4 )
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* obp = &ptr[ cnt * 4 + 3 ];

                *obp = (uchar)( (float)*obp * val );
            }
        }
        else
        {
            #pragma omp parallel for
            for( OMPSIZE_T cnt=0; cnt<img_sz; cnt++ )
            {
                uchar* obp = &ptr[ cnt * 3 ];

                for( OMPSIZE_T rpt=0; rpt<3; rpt++ )
                {
                    obp[rpt] = (uchar)( (float)obp[rpt] * val );
                }
            }
        }
        
        return true;
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

// Now it handles alpha channels !
Fl_RGB_Image* fl_imgtk::makegradation_h( unsigned w, unsigned h, ulong col1, ulong col2, bool dither )
{
    // sens alpha channel using.
    unsigned d = 4;
    uchar alpha1 = col1 & 0x000000FF;
    uchar alpha2 = col2 & 0x000000FF;

    if ( ( alpha1 == 0 ) || ( alpha2 == 0 ) )
    {
        d = 3;
        
        alpha1 = 0xFF;
        alpha2 = 0xFF;
    }

    uchar* buffer = new unsigned char[ w * h * d ];
    if ( buffer != NULL )
    {
        float downscale_f = 2559.0f / ( float( h ) * 10.0f );
        uchar ref_c1r   = ( col1 & 0xFF000000 ) >> 24;
        uchar ref_c1g   = ( col1 & 0x00FF0000 ) >> 16;
        uchar ref_c1b   = ( col1 & 0x0000FF00 ) >> 8;
        uchar ref_c1a   = alpha1;
        uchar ref_c2r   = ( col2 & 0xFF000000 ) >> 24;
        uchar ref_c2g   = ( col2 & 0x00FF0000 ) >> 16;
        uchar ref_c2b   = ( col2 & 0x0000FF00 ) >> 8;
        uchar ref_c2a   = alpha2;

        #pragma omp parallel for
        for ( long cy=0; cy<(long)h; cy++ )
        {
            float col1g     = 0.0f;
            float col2g     = 0.0f;
            float randfv    = 0.0f;
            uchar fill_r    = 0;
            uchar fill_g    = 0;
            uchar fill_b    = 0;
            uchar fill_a    = 0;
            uchar* sbuf = &buffer[ cy * w * d ];

            col1g = ( ( float( h*10 ) - float( cy*10 )  ) * downscale_f ) / 2559.0f;
            col2g = 1.0f - col1g;

            fill_r = ( (float)ref_c1r * col1g ) + ( (float)ref_c2r * col2g );
            fill_g = ( (float)ref_c1g * col1g ) + ( (float)ref_c2g * col2g );
            fill_b = ( (float)ref_c1b * col1g ) + ( (float)ref_c2b * col2g );
            if ( d > 3 )
            {
                fill_a = ( (float)ref_c1a * col1g ) + ( (float)ref_c2a * col2g );
            }

            for ( long cx=0; cx<(long)w; cx++ )
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
                    if ( d > 3 )
                    {
                        sbuf[3] = fill_a - uchar( ref_c2a * col2g * randfv );
                    }
                }
                else
                {
                    sbuf[0] = fill_r;
                    sbuf[1] = fill_g;
                    sbuf[2] = fill_b;
                    if ( d > 3 )
                    {
                        sbuf[3] = fill_a;
                    }
                }
                sbuf += d;
            }
        }
#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buffer, w, h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buffer, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
    }

    return NULL;
}

Fl_RGB_Image* fl_imgtk::makegradation_v( unsigned w, unsigned h, ulong col1, ulong col2, bool dither )
{
    // sens alpha channel using.
    unsigned d = 4;
    uchar alpha1 = col1 & 0x000000FF;
    uchar alpha2 = col2 & 0x000000FF;

    if ( ( alpha1 == 0 ) || ( alpha2 == 0 ) )
    {
        d = 3;
        alpha1 = 0xFF;
        alpha2 = 0xFF;
    }

    uchar* buffer = new unsigned char[ w * h * d ];
    if ( buffer != NULL )
    {
        float downscale_f = 2559.0f / ( float( w ) * 10.0f );
        uchar ref_c1r   = ( col1 & 0xFF000000 ) >> 24;
        uchar ref_c1g   = ( col1 & 0x00FF0000 ) >> 16;
        uchar ref_c1b   = ( col1 & 0x0000FF00 ) >> 8;
        uchar ref_c1a   = alpha1;
        uchar ref_c2r   = ( col2 & 0xFF000000 ) >> 24;
        uchar ref_c2g   = ( col2 & 0x00FF0000 ) >> 16;
        uchar ref_c2b   = ( col2 & 0x0000FF00 ) >> 8;
        uchar ref_c2a   = alpha2;

        OMPSIZE_T cx = 0;
        OMPSIZE_T cy = 0;
        
        #pragma omp parallel for private( cy )
        for ( cx=0; cx<w; cx++ )
        {
            float col1g     = 0.0f;
            float col2g     = 0.0f;
            float randfv    = 0.0f;
            uchar fill_r    = 0;
            uchar fill_g    = 0;
            uchar fill_b    = 0;
            uchar fill_a    = 0;
            
            col1g = ( ( float( w*10 ) - float( cx*10 )  ) * downscale_f ) / 2559.0f;
            col2g = 1.0f - col1g;

            fill_r = ( (float)ref_c1r * col1g ) + ( (float)ref_c2r * col2g );
            fill_g = ( (float)ref_c1g * col1g ) + ( (float)ref_c2g * col2g );
            fill_b = ( (float)ref_c1b * col1g ) + ( (float)ref_c2b * col2g );
            if ( d > 3 )
            {
                fill_a = ( (float)ref_c1a * col1g ) + ( (float)ref_c2a * col2g );
            }

            for ( cy=0; cy<h; cy++ )
            {
                uchar* sbuf = &buffer[ ( ( w * cy ) + cx ) * d ];

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
                    if ( d > 3 )
                    {
                        sbuf[3] = fill_a - uchar( ref_c2a * col2g * randfv );
                    }
                }
                else
                {
                    sbuf[0] = fill_r;
                    sbuf[1] = fill_g;
                    sbuf[2] = fill_b;
                    if ( d > 3 )
                    {
                        sbuf[3] = fill_a;
                    }
                }
                sbuf += d;
            }
        }

#if defined(FLIMGTK_IMGBUFF_OWNALLOC)
        Fl_RGB_Image* newimg = new Fl_RGB_Image( buffer, w, h, d );
        if ( newimg != NULL )
        {
            newimg->alloc_array = 1;
            return newimg;
        }
#else
        return new Fl_RGB_Image( buffer, w, h, d );
#endif /// of #if defined(FLIMGTK_IMGBUFF_OWNALLOC)
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

inline void fl_imgtk_dla_plot( Fl_RGB_Image* img, int x, int y, ulong col, float br )
{
    float alpha = (float)( col & 0x000000FF ) / 255.f;
    
    if ( ( br <= 0.f ) || ( alpha == 0.f ) )
        return;

    uchar col_r = ( col & 0xFF000000 ) >> 24;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_b = ( col & 0x0000FF00 ) >> 8;

    int w = img->w();
    int h = img->h();
    int d = img->d();

    if ( d < 3 )
    {
        unsigned col_avr = ( col_r + col_g + col_b ) / 3;
        col_r = col_g = col_b = (uchar)(col_avr & 0x000000FF );
    }

    br *= alpha;

    if ( ( ( x >= 0 ) && ( y >= 0 ) )
        && ( ( x < w ) && ( y < h ) ) )
    {
        unsigned pos  = ( y * w + x ) * d;
        float revbr = 1.0f - br;

        uchar* ptrimg = (uchar*)img->data()[0];

        if ( d >= 3 ) /// up to RGB and RGBA
        {
            ptrimg[ pos + 0 ] = ( ptrimg[ pos + 0 ] * revbr ) + ( (float)col_r * br );
            ptrimg[ pos + 1 ] = ( ptrimg[ pos + 1 ] * revbr ) + ( (float)col_g * br );
            ptrimg[ pos + 2 ] = ( ptrimg[ pos + 2 ] * revbr ) + ( (float)col_b * br );
            if ( d > 3 )
            {
                ptrimg[ pos + 3 ] = (uchar)( ( ptrimg[ pos + 3 ] * revbr ) + ( alpha * 255.f ) * br );
            }
        }
        else /// up to Y or Y+A
        {
            ptrimg[ pos ] = ( ptrimg[ pos ] * revbr ) + ( (float)col_r * br );
            if ( d > 1 )
            {
                ptrimg[ pos + 1 ] = (uchar)( ( ptrimg[ pos + 1 ] * revbr ) + ( alpha * 255.f ) * br );
            }
        }
    }
}

#define fl_imgtk_ipart(X)   ((int)(X))
#define fl_imgtk_roundi(X)  ((int)(((float)(X))+0.5f))
#define fl_imgtk_fpart(X)   (((float)(X))-(float)fl_imgtk_ipart(X))
#define fl_imgtk_rfpart(X)  (1.0f - fl_imgtk_fpart(X))
#ifdef _MSC_VER /// godam M$VC !!
    #define fl_imgtk_swap(a, b) { auto tmp = a; a = b; b = tmp; }
#else
    #define fl_imgtk_swap(a, b) { __typeof__(a) tmp;  tmp = a; a = b; b = tmp; }
#endif

void fl_imgtk::draw_smooth_line( Fl_RGB_Image* img, int x1, int y1, int x2, int y2, ulong col )
{
    float dx = (float)x2 - (float)x1;
    float dy = (float)y2 - (float)y1;

    if ( fabs(dx) > fabs(dy) )
    {
        if ( x2 < x1 )
        {
            fl_imgtk_swap( x1, x2 );
            fl_imgtk_swap( y1, y2 );
        }

        float gradient = dy / dx;
        float xend = fl_imgtk_roundi( x1 );
        float yend = (float)y1 + gradient*( xend - x1 );
        float xgap = fl_imgtk_rfpart( xend );

        int xpxl1 = xend;
        int ypxl1 = fl_imgtk_ipart( yend );

        fl_imgtk_dla_plot( img, xpxl1, ypxl1, col, 
                           fl_imgtk_rfpart(yend)*xgap );
        fl_imgtk_dla_plot( img, xpxl1, ypxl1+1, col, 
                           fl_imgtk_fpart(yend)*xgap );

        float intery = yend + gradient;

        xend = fl_imgtk_roundi( x2 );
        yend = y2 + gradient*( xend - x2 );
        xgap = fl_imgtk_fpart( xend );

        int xpxl2 = xend;
        int ypxl2 = fl_imgtk_ipart( yend );

        fl_imgtk_dla_plot( img, xpxl2, ypxl2, col, 
                           fl_imgtk_rfpart(yend) * xgap );
        fl_imgtk_dla_plot( img, xpxl2 + 1, ypxl2, col, 
                           fl_imgtk_fpart(yend) * xgap );

        for( int x=xpxl1+1; x<=xpxl2; x++ )
        {
            fl_imgtk_dla_plot( img, x, fl_imgtk_ipart(intery), col, 
                               fl_imgtk_rfpart(intery) );
            fl_imgtk_dla_plot( img, x, fl_imgtk_ipart(intery) + 1, col, 
                               fl_imgtk_fpart(intery) );
            intery += gradient;
        }
    }
    else
    {
        if ( y2 < y1 )
        {
            fl_imgtk_swap( x1, x2 );
            fl_imgtk_swap( y1, y2 );
        }

        float gradient = dx / dy;
        float yend = fl_imgtk_roundi( y1 );
        float xend = (float)x1 + gradient*( yend - y1 );
        float ygap = fl_imgtk_rfpart( yend );

        int ypxl1 = yend;
        int xpxl1 = fl_imgtk_ipart(xend);

        fl_imgtk_dla_plot( img, xpxl1, ypxl1, col, 
                           fl_imgtk_rfpart(xend)*ygap );
        fl_imgtk_dla_plot( img, xpxl1, ypxl1+1, col, 
                           fl_imgtk_fpart(xend)*ygap );

        float interx = xend + gradient;

        yend = fl_imgtk_roundi( y2 );
        xend = x2 + gradient*( yend - y2 );
        ygap = fl_imgtk_fpart( yend );

        int ypxl2 = yend;
        int xpxl2 = fl_imgtk_ipart( xend );

        fl_imgtk_dla_plot( img, xpxl2, ypxl2, col, fl_imgtk_rfpart(xend) * ygap );
        fl_imgtk_dla_plot( img, xpxl2, ypxl2 + 1, col, fl_imgtk_fpart(xend) * ygap );

        for( int y=ypxl1+1; y<=ypxl2; y++ )
        {
            fl_imgtk_dla_plot( img, fl_imgtk_ipart(interx), y, col, fl_imgtk_rfpart(interx) );
            fl_imgtk_dla_plot( img, fl_imgtk_ipart(interx) + 1, y, col, fl_imgtk_fpart(interx) );
            interx += gradient;
        }
    }
}

inline float _capsule_sdf(float px, float py, float ax, float ay, float bx, float by, float r) 
{
    float pax = px - ax, pay = py - ay, bax = bx - ax, bay = by - ay;
    float h   = fmaxf(fminf((pax * bax + pay * bay) / (bax * bax + bay * bay), 1.0f), 0.0f);
    float dx  = pax - bax * h, dy = pay - bay * h;
    return sqrtf(dx * dx + dy * dy) - r;
}

// reference : https://github.com/miloyip/line
// algorithm replaced.
void fl_imgtk::draw_smooth_line_ex( Fl_RGB_Image* img, int x1, int y1, int x2, int y2, float wd, ulong col )
{
    if ( img == NULL )
        return;
    
    float r = wd / 2;
    
    if ( r < 0.1f )
        r = 0.1f;
        
    long _x0 = (long)floorf(fminf(x1, x2) - r);
    long _x1 = (long) ceilf(fmaxf(x1, x2) + r);
    long _y0 = (long)floorf(fminf(y1, y2) - r);
    long _y1 = (long) ceilf(fmaxf(y1, y2) + r);

    #pragma omp parallel for
    for (long y = _y0; y <= _y1; y++)
    {
        for (long x = _x0; x <= _x1; x++)
        {
            fl_imgtk_dla_plot( img, x, y, col, 
                               fmaxf(fminf(0.5f - _capsule_sdf(x, y, x1, y1, x2, y2, r), 1.0f), 0.0f) );
        }
    }
}

#if defined(FAST_RGBA_IGNORE_ALPHA)
#define fl_imgtk_putpixel( _buff_,_x_,_w_,_y_,_d_,_r_,_g_,_b_,_a_ ) \
        uchar* _putptr_ = &_buff_[ ( ( _y_ * _w_ ) + _x_ ) * _d_ ];\
        if((uchar)_a_==0xFF)\
        {_putptr_[0]=_r_; _putptr_[1]=_g_; _putptr_[2]=_b_;}\
        else\
        {float _ar_=(float)_a_/255.f; float _rar_=1.f-_ar_;\
         _putptr_[0]=(uchar)(((float)_putptr_[0]*_rar_) + ((float)_r_*_ar_));\
         _putptr_[1]=(uchar)(((float)_putptr_[1]*_rar_) + ((float)_g_*_ar_));\
         _putptr_[2]=(uchar)(((float)_putptr_[2]*_rar_) + ((float)_b_*_ar_));}
#else
inline void fl_imgtk_putpixel( uchar* buff, \
                               unsigned x, unsigned w, unsigned y, unsigned d, \
                               unsigned r, unsigned g, unsigned b, unsigned a )
{
    if ( a == 0 )
        return;

    uchar* putptr = &buff[ ( ( y * w ) + x ) * d ];
    
    float ar  = (float)a / 255.f; 
    float rar = 1.f - ar;

    if ( d >= 3 )
    {
        putptr[0]=(uchar)(((float)putptr[0]*rar) + ((float)r*ar));
        putptr[1]=(uchar)(((float)putptr[1]*rar) + ((float)g*ar));
        putptr[2]=(uchar)(((float)putptr[2]*rar) + ((float)b*ar));
        if ( d > 3 )
        {
            putptr[3]=(uchar)(((float)putptr[3]*rar)+(ar*255.f)*ar);
        }
    }
    else
    {
        float colavr = (((float)r*ar) + ((float)g*ar) + ((float)b*ar))/3.f;
        putptr[0] = (uchar)(((float)putptr[0]*rar) + colavr );
        if ( d > 1 )
        {
            putptr[1]=(uchar)(((float)putptr[1]*rar)+(ar*255.f)*ar);
        }
    }
}
#endif // of defined(FAST_RGBA_IGNORE_ALPHA)

void fl_imgtk::draw_line( Fl_RGB_Image* img, int x1, int y1, int x2, int y2, ulong col )
{
    if ( img == NULL )
    {
        return;
    }
    
    uchar* putbuff = (uchar*)img->data()[0];
    uchar col_b = ( col & 0x0000FF00 ) >> 8;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_r = ( col & 0xFF000000 ) >> 24;
    uchar col_a = ( col & 0x000000FF );

    if ( col_a == 0 )
        return;

    uchar img_d = img->d();
    int   img_w = img->w();
    int   img_h = img->h();

    int _x0   = x1;
    int _x1   = x2;
    int _y0   = y1;
    int _y1   = y2;

    int inc1  = 0;
    int inc2  = 0;
    int cnt   = 0;
    int y_adj = 0;
    int dy    = 0;
    int dx    = 0;
    int x_adj = 0;

    if ( _x0 == _x1 )
    {
        if ( _y0 > _y1 )
        {
            swap( _y0, _y1 );
        }

        int cnt = _y1 - _y0 + 1;

        while( cnt-- )
        {
            if ( ( ( _x0 >= 0 ) && ( _x0 < img_w ) ) && 
                 ( ( _y0 + cnt >= 0 ) && ( _y0 + cnt < img_h ) ) )
            {
                fl_imgtk_putpixel( putbuff, _x0, img_w, (_y0 + cnt), 
                                   img_d, col_r, col_g, col_b, col_a );
            }
        }
    }
    else
    {
        if ( _y0 == _y1 )
        {
            if ( _x0 > _x1 )
            {
                swap( _x0, _x1 );
            }

            dx = _x1 - _x0 + 1;

            for( int cnt=0; cnt<dx; cnt++ )
            {
                if ( ( ( _x0 + cnt >= 0 ) && ( _x0 + cnt < img_w ) ) && 
                     ( ( _y0 >= 0 ) && ( _y0 < img_h ) ) )
                {
                    fl_imgtk_putpixel( putbuff, (_x0 + cnt), img_w, _y0,
                                       img_d, col_r, col_g, col_b, col_a );
                }
            }
        }
        else
        {
            dy = _y1 - _y0;
            dx = _x1 - _x0;

            if ( abs( dy ) < abs( dx ) )
            {
                if ( _x0 > _x1 )
                {
                    swap( _x0, _x1 );
                    swap( _y0, _y1 );
                }

                _y0--;
                _y1--;
                _x0--;
                _x1--;

                dy = _y1 - _y0;
                dx = _x1 - _x0;

                if ( dy < 0 )
                {
                    dy    = -dy;
                    y_adj = -1;
                }
                else
                if ( dy > 0 )
                {
                    y_adj = 1;
                }
                else
                {
                    y_adj = 0;
                }

                inc1 = dy << 1;
                inc2 = ( dy - dx ) << 1;
                cnt  = ( dy << 1 ) - dx;

                dx++;
                int py = y_adj;

                while ( dx-- )
                {
                    if ( ( ( _x0 + dx >= 0 ) && ( _x0 + dx < img_w ) ) && 
                         ( ( _y1 - py >= 0 ) && ( _y1 - py < img_h ) ) )
                    {
                        fl_imgtk_putpixel( putbuff, 
                                           _x0 + dx, 
                                           img_w, 
                                           (_y1 - py), 
                                           img_d, 
                                           col_r, col_g, col_b, col_a );
                    }

                    if ( cnt >= 0 )
                    {
                        cnt += inc2;
                        py  += y_adj;
                    }
                    else
                    {
                        cnt += inc1;
                    }
                }
            }
            else
            {
                if ( _y0 > _y1 )
                {
                    swap( _x0, _x1 );
                    swap( _y0, _y1 );
                }

                _x0--;
                _y0--;
                _x1--;
                _y1--;

                dy = _y1 - _y0;
                dx = _x1 - _x0;

                if ( dx < 0)
                {
                    dx    = -dx;
                    x_adj = -1;
                }
                else
                if ( dx > 0 )
                {
                    x_adj = 1;
                }
                else
                {
                    x_adj = 0;
                }

                inc1 = dx << 1;
                inc2 = ( dx - dy ) << 1;
                cnt  = ( dx << 1 ) - dy;

                dy++;
                int px = x_adj;

                while ( dy-- )
                {
                    if ( ( ( _x0 + px >= 0 ) && ( _x0 + px < img_w ) ) && 
                         ( ( _y1 - dy >= 0 ) && ( _y1 - dy < img_h ) ) )
                    {                   
                        fl_imgtk_putpixel( putbuff, 
                                           _x0 + px, 
                                           img_w, 
                                           (_y1 - dy), 
                                           img_d, 
                                           col_r, col_g, col_b, col_a );
                    }

                    if ( cnt >= 0 )
                    {
                        cnt += inc2;
                        px  += x_adj;
                    }
                    else
                    {
                        cnt += inc1;
                    }
                }
            }
        }
    }
}

void fl_imgtk::draw_rect( Fl_RGB_Image* img, int x, int y, int w, int h, ulong col )
{
    if ( ( w < 0 ) || ( h < 0 ) ) return;

    if ( x < 0 )
    {
        w += x;
        x = 0;
    }

    if ( y < 0 )
    {
        h += y;
        y = 0;
    }

    int x1 = x;
    int y1 = y;
    int x2 = x + w;
    int y2 = y + h;

    if ( ( x2 < 1 ) && ( y2 < 1 ) )
        return;
          
    // =
    draw_line( img, x1, y1, x2, y1, col );
    draw_line( img, x1, y2, x2, y2, col );

    if ( ( x2 - x1 ) > 1 )
    {
        y1++;
        if ( ( x2 - x1 ) > 2 )
        {
            y2--;
        }
    }

    // ||
    draw_line( img, x1, y1, x1, y2, col );
    draw_line( img, x2, y1, x2, y2, col );
}

void fl_imgtk::draw_fillrect( Fl_RGB_Image* img, int x, int y, int w, int h, ulong col )
{
    if ( ( w < 0 ) || ( h < 0 ) ) return;

    if ( x < 0 )
    {
        w += x;
        x = 0;
    }

    if ( y < 0 )
    {
        h += y;
        y = 0;
    }

    uchar* putbuff = (uchar*)img->data()[0];
    uchar col_b = ( col & 0x0000FF00 ) >> 8;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_r = ( col & 0xFF000000 ) >> 24;
    uchar col_a = ( col & 0x000000FF );
    
    OMPSIZE_T img_d = img->d();
    OMPSIZE_T img_w = img->w();
    OMPSIZE_T img_h = img->h();

    OMPSIZE_T cym = y+h;
    OMPSIZE_T cxm = x+w;

    if ( cxm > img_w )
        cxm = img_w - x;

    if ( cym > img_h )
        cym = img_h - y;

    #pragma omp parallel for
    for( OMPSIZE_T cy=y; cy<cym; cy++ )
    {
        for( OMPSIZE_T cx=x; cx<cxm; cx++ )
        {
            fl_imgtk_putpixel( putbuff,
                               cx, img_w, cy, img_d,
                               col_r, col_g, col_b, col_a );
        }
    }
}

bool fl_imgtk_sortcondition (int i,int j)
{
    return ( i < j );
}

void fl_imgtk::draw_polygon( Fl_RGB_Image* img, const fl_imgtk::vecpoint* points, unsigned pointscnt, ulong col )
{
    if ( img == NULL )
        return;
    
    if ( ( points == NULL ) || ( pointscnt < 3 ) )
        return;
        
    uchar col_r = ( col & 0xFF000000 ) >> 24;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_b = ( col & 0x0000FF00 ) >> 8;
    uchar col_a = ( col & 0x000000FF );
    
    uchar* ptrimg = (uchar*)img->data()[0];
    unsigned img_w = img->w();
    unsigned img_h = img->h();
    unsigned img_d = img->d();
    unsigned max_y = img_h;
    unsigned max_x = img_w;

    vector< double > node_x;

    for( unsigned cur_y=0; cur_y<max_y; cur_y++ )
    {
        unsigned    ptsz        = pointscnt;
        unsigned    rcnt        = ptsz - 1;
        unsigned    node_count  = 0;

        for( unsigned cnt=0; cnt<ptsz; cnt++ )
        {
            double pt_x  = (double)points[ cnt ].x;
            double pt_y  = (double)points[ cnt ].y;
            double pt_rx = (double)points[ rcnt ].x;
            double pt_ry = (double)points[ rcnt ].y;
            double dc_y  = (double)cur_y;

            if ( ( ( pt_y  < dc_y ) && ( pt_ry >= dc_y ) ) ||
                 ( ( pt_ry < dc_y ) && ( pt_y  >= dc_y ) ) )
            {
                double newv = pt_x + ( dc_y - pt_y ) /
                              ( pt_ry - pt_y ) * ( pt_rx - pt_x );

                node_x.push_back( newv );
            }

            rcnt = cnt;
        }

        OMPSIZE_T node_x_sz = (OMPSIZE_T)node_x.size();

        // sort nodes ..
        if ( node_x_sz > 1 )
        {
            sort( node_x.begin(),
                  node_x.begin() + node_x_sz,
                  fl_imgtk_sortcondition );
#if (DEBUG_DRAWPOLYGON==1)
            for( unsigned ncnt=0; ncnt<node_x_sz; ncnt++ )
            {
                printf("%.2f, ", node_x[ncnt] );
            }
            printf("\n");
#endif /// of DEBUG_DRAWPOLYGON
        }

        #pragma parallel for
        for( OMPSIZE_T dcnt=0; dcnt<node_x_sz; dcnt+=2 )
        {
            if ( node_x[dcnt] >= max_x )
                break;

            if ( node_x[dcnt+1] > 0 )
            {
                if ( node_x[dcnt] < 0 )
                {
                    node_x[dcnt] = 0;
                }

                if ( node_x[dcnt+1] > max_x )
                {
                    node_x[dcnt+1] = max_x;
                }

                for( OMPSIZE_T cur_x=node_x[dcnt]; cur_x<=node_x[dcnt+1]; cur_x++ )
                {
                    fl_imgtk_putpixel( ptrimg,
                                       cur_x, img_w, cur_y, img_d,
                                       col_r, col_g, col_b, col_a );
                }
            }
        }

        node_x.clear();

    } /// of for( y .. )

    img->uncache();
}

void fl_imgtk::draw_2xaa_polygon( Fl_RGB_Image* img, const vecpoint* points, unsigned pointscnt, ulong col )
{
    if ( img == NULL )
        return;
    
    if ( ( points == NULL ) || ( pointscnt < 3 ) )
        return;
    
    if ( img->d() < 3 )
        return;
    
    uchar col_r = ( col & 0xFF000000 ) >> 24;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_b = ( col & 0x0000FF00 ) >> 8;
    uchar col_a = ( col & 0x000000FF );
    
    // get max coordinates ...
    unsigned max_x = 0;
    unsigned max_y = 0;
    int      min_x = 0xFFFFFFFF;
    int      min_y = 0xFFFFFFFF;
    
    for ( size_t pcnt=0; pcnt<pointscnt; pcnt++ )
    {
        if ( points[pcnt].x > max_x ) max_x = (unsigned)points[pcnt].x;
        if ( points[pcnt].y > max_y ) max_y = (unsigned)points[pcnt].y;
        if ( points[pcnt].x < min_x ) min_x = (unsigned)points[pcnt].x;
        if ( points[pcnt].y < min_y ) min_y = (unsigned)points[pcnt].y;
    }

    if ( min_x < 0 )
    {
        max_x += (unsigned)(-min_x);
    }
    
    if ( min_y < 0 )
    {
        max_y += (unsigned)(-min_y);
    }

    // this is a trick, over sized drawn to downscale to 100%.
    float mulrat_w = 1.8f;
    float mulrat_h = 1.8f;
    float divrat_w = 1.f / mulrat_w;
    float divrat_h = 1.f / mulrat_h;

    if ( ( max_x > 0 ) && ( max_y > 0 ) )
    {
        unsigned img_w = max_x * mulrat_w;
        unsigned img_h = max_y * mulrat_h;
        
        vecpoint* dblpts = new vecpoint[ pointscnt ];
        
        if ( dblpts == NULL )
            return;
            
        for ( size_t cnt=0; cnt<pointscnt; cnt++ )
        {
            dblpts[cnt].x = points[cnt].x * mulrat_w;
            dblpts[cnt].y = points[cnt].y * mulrat_h;
        }
        
        Fl_RGB_Image* imgdbl = makeanempty( img_w, img_h, 4, 0x00000000 );
        
        if ( imgdbl == NULL )
        {
            delete[] dblpts;
            return;
        }
        
        draw_polygon( imgdbl, dblpts, pointscnt, col );

        delete[] dblpts;

        // downscale to make antialiased --
        Fl_RGB_Image* imgdsl = rescale( imgdbl, 
                                        img_w*divrat_w, img_h*divrat_h, 
                                        BILINEAR );
        
        if ( imgdsl == NULL )
        {
            discard_user_rgb_image( imgdbl );
            return;
        }
        
        drawonimage( img, imgdsl, min_x, min_y );
        
        discard_user_rgb_image( imgdbl );
        discard_user_rgb_image( imgdsl );
    }
    
    img->uncache();
}

void fl_imgtk::discard_user_rgb_image( Fl_RGB_Image* &img )
{
    if( img != NULL )
    {
#if !defined(FLIMGTK_IMGBUFF_OWNALLOC)
        if ( ( img->array != NULL ) && ( img->alloc_array == 0 ) )
        {
            delete[] img->array;
        }
#endif
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
