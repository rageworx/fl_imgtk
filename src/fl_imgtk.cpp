#ifdef _MSC_VER
#pragma warning(disable : 4018)  
#pragma warning(disable : 4068)  
#pragma warning(disable : 4244)  
#pragma warning(disable : 4996)  
#endif

#include "fl_imgtk.h"
#include "fl_imgtk_minmax.h"

#include <FL/Fl_Widget.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Image_Surface.H>
#include <FL/fl_draw.H>

#include "fl_smimg.h"

#include "stdint.h"

#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

#ifdef USING_OMP
#include <omp.h>
#endif /// of USING_OMP

////////////////////////////////////////////////////////////////////////////////

#define FLOAT_PI        3.141592654f
#define FLOAT_PI2X      6.28318530718f

#define FLIMGTK_BI_RGB       0  /// No compression - straight BGR data
#define FLIMGTK_BI_RLE8      1  /// 8-bit run-length compression
#define FLIMGTK_BI_RLE4      2  /// 4-bit run-length compression
#define FLIMGTK_BI_BITFIELDS 3  /// RGB bitmap with RGB masks

#define fl_imgtk_degree2f( _x_ )      ( ( _x_ / 360.f ) * FLOAT_PI2X )
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


Fl_RGB_Image* fl_imgtk::makeanempty( unsigned w, unsigned h, unsigned d, ulong color )
{
    if ( ( w > 0 ) && ( h > 0 ) && ( d >= 3 ) )
    {
        ulong  resz   = w * h;
        ulong  datasz = resz * d;
        uchar* pdata  = new uchar[ datasz ];
        
        uchar ref_r   = ( color & 0xFF000000 ) >> 24;
        uchar ref_g   = ( color & 0x00FF0000 ) >> 16;
        uchar ref_b   = ( color & 0x0000FF00 ) >> 8;
        uchar ref_a   = ( color & 0x000000FF );

        uchar carray[4] = { ref_r, ref_g, ref_b, ref_a };
        
        if ( pdata != NULL )
        {
            #pragma omp parallel for
            for( ulong cnt=0; cnt<resz; cnt++ )
            {
                memcpy( &pdata[ cnt * d ], &carray[0], d );
            }
            
            return new Fl_RGB_Image( pdata, w, h, d );
        }
    }
    
    return NULL;
}

uint16_t flimgtk_memread_word( const char* buffer, unsigned* que = NULL )
{
    uint16_t retus = 0;
    memcpy( &retus, buffer, 2 );
    
    if ( que != NULL )
    {
        *que += 2;
    }
    
    return retus;
}

uint32_t flimgtk_memread_dword( const char* buffer, unsigned* que = NULL )
{
    uint32_t retui = 0;
    memcpy( &retui, buffer, 4 );
    
    if ( que != NULL )
    {
        *que += 4;
    }
    
    return retui;
}

int flimgtk_memread_int( const char* buffer, unsigned* que = NULL )
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
        buffque += 4;   /// skip size.
        buffque += 2;   /// skip reserved something.
        buffque += 2;
        
        int info_size;      /// Size of info header
        uint16_t depth;     /// Depth of image (bits)
        int bDepth = 3;     /// Depth of image (bytes)
        uint16_t compression; /// Type of compression
        uint32_t colors_used; /// Number of colors used
        int x, y;           /// Looping vars
        int32_t color = 0;  /// Color of RLE pixel
        int repcount;       /// Number of times to repeat
        int temp;           /// Temp. Color or Index
        int align;          /// Alignment bytes
        uint32_t dataSize;    /// number of bytes in image data set
        int row_order=-1;   /// 1 = normal;  -1 = flipped row order
        int start_y;        /// Beginning Y
        int end_y;          /// Ending Y
        int offbits;        // Offset to image data
        uchar bit;          /// Bits in image
        uchar byte;         /// Bytes in image
        uchar*ptr;          /// Pointer into pixels
        uchar colormap[256][3];/// Colormap
        uchar havemask = 0; /// Single bit mask follows image data
        int use_5_6_5 = 0;  /// Use 5:6:5 for R:G:B channels in 16 bit images
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
        
        return new Fl_RGB_Image( array, w, h, bDepth );
    }
        
    return NULL;
}

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

        //#pragma omp parallel for private( cnth )
        for( cntw=new_w; cntw-- != 0; )
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

        //#pragma omp parallel for private( cnth )
        for( cntw=0; cntw<new_w; cntw++ )
        {
            for( cnth=new_h; cnth-- != 0; )
            {
                unsigned pos1 = ( new_w * cnth + cntw ) * d;
                unsigned pos2 = ( src_w * cntw + new_h - cnth - 1 ) * d;

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
// http://www.codeguru.com/cpp/g-m/gdi/article.php/c3693/ ...
// Rotate-a-Bitmap-at-Any-Angle-Without-GetPixelSetPixel.htm
// Maybe, it will slowly works if not enable openMP.
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

    memset( obuff, 0, dstW * dstH * img_d );

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

            if ( ( (orgX >= 0) && (orgY >= 0) ) && \
                 ( (orgX < img_w-1) && (orgY < img_h-1) ) )
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
            return false;

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
            return false;

        if ( ( kfc->w == 0 ) || ( kfc->h == 0 ) || ( kfc->msz == 0 ) || ( kfc->m == NULL ) )
            return false;

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

bool fl_imgtk_gen_lowfreq( uchar** out, uchar* src, unsigned w, unsigned h, unsigned d, unsigned sz )
{
    if ( src == NULL )
        return false;

    if ( ( w == 0 ) || ( h == 0 ) || ( d < 3 ) )
        return false;

    if ( sz < 2 )
        return false;

    long startpos = sz / 2;
    long endposw  = w - startpos;
    long endposh  = h - startpos;

    uchar* outbuff = new uchar[ w * h * d ];

    if ( outbuff == NULL )
        return false;

    long cnty;
    long cntx;

    #pragma omp parallel for private( cntx )
    for( cnty=startpos; cnty<endposh; cnty++ )
    {
        for( cntx=startpos; cntx<endposw; cntx++ )
        {
            unsigned sum[4] = {0};

            for( long ry = -startpos; ry<(startpos+1); ry++ )
            {
                for( long rx = -startpos; rx<(startpos+1); rx++ )
                {
                    for( long rpt=0; rpt<d; rpt++ )
                    {
                        sum[rpt] += src[ ( ( cnty + ry ) * w + ( cntx + rx ) ) * d + rpt ];
                    }
                }
            }

            for( long rpt=0; rpt<d; rpt++ )
            {
                outbuff[ ( cnty * w + cntx ) * d + rpt ] \
                 = MIN( 255, sum[rpt] / ( sz * sz ) );
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
        unsigned imgsz = img->w() * img->h() * img->d();
        uchar* outbuff = new uchar[ imgsz ];

        if ( outbuff == NULL )
            return NULL;

        uchar* rbuff  = (uchar*)img->data()[0];
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

        unsigned cnty;
        unsigned cntx;
        unsigned cntw = img->w();
        unsigned cnth = img->h();
        unsigned mgnx = margin;
        unsigned mgny = margin;
        unsigned mgnw = img->w() - ( margin * 2 );
        unsigned mgnh = img->h() - ( margin * 2 );

        #pragma omp parallel for private( cntx )
        for( cnty=0; cnty<cnth; cnty++ )
        {
            for( cntx=0; cntx<cntw; cntx++ )
            {
                for( unsigned rpt=0; rpt<img->d(); rpt++ )
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
            }
        }

        delete[] lfimg5;
        delete[] lfimg9;

        return new Fl_RGB_Image( outbuff, img->w(), img->h(), img->d() );
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

        unsigned cnty;
        unsigned cntx;
        unsigned mgnx = margin;
        unsigned mgny = margin;
        unsigned mgnw = img->w() - ( margin * 2 );
        unsigned mgnh = img->h() - ( margin * 2 );

        #pragma omp parallel for private( cntx )
        for( cnty=mgny; cnty<mgnh; cnty++ )
        {
            for( cntx=mgnx; cntx<mgnw; cntx++ )
            {
                for( unsigned rpt=0; rpt<img->d(); rpt++ )
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

            Fl_RGB_Image* retimg =  new Fl_RGB_Image( widgetbuff,
                                                      cwin_w,
                                                      cwin_h,
                                                      3 );
                                                      
            if ( retimg == NULL )
            {
                delete[] widgetbuff;
            }
            
            return retimg;
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

        if ( ( rbuff != NULL ) && ( obuff != NULL ) )
        {
            unsigned srcw = rsx + w;
            unsigned srch = rsy + h;
            unsigned cnty;
            unsigned cntx;
            unsigned putx;
            unsigned puty;

            #pragma omp parellel for private( cntx )
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
#ifdef DEBUG
            Fl_RGB_Image* retimg = new Fl_RGB_Image( obuff, rw, rh , sd );
            return retimg;
#else
            return new Fl_RGB_Image( obuff, rw, rh , sd );
#endif      
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

                float a_opa = MIN( 1.f, alpha );
                
                float w_r = (float)wptr[0] / 255.f;
                float w_g = (float)wptr[1] / 255.f;
                float w_b = (float)wptr[2] / 255.f;
                float w_a = 1.0f;
                
                float r_r = (float)rptr[0] / 255.f;
                float r_g = (float)rptr[1] / 255.f;
                float r_b = (float)rptr[2] / 255.f;
                float r_a = 1.0f;

                if ( bd == 4 )
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
                    
                float rTemp = r_r * r_a + w_r * w_a * ( 1.0f - r_a );
                float gTemp = r_g * r_a + w_g * w_a * ( 1.0f - r_a );
                float bTemp = r_b * r_a + w_b * w_a * ( 1.0f - r_a );

                r_a = w_a + ( 1.0f - w_a ) * r_a;

                if ( r_a == 0.0f )
                {
                    if ( a_opa > 0.0f )
                    {
                        rTemp = r_r;
                        gTemp = r_g;
                        bTemp = r_b;
                    }
                    else
                    {
                        rTemp = w_r;
                        gTemp = w_g;
                        bTemp = w_b;
                    }
                }
                else
                {
                    rTemp /= r_a;
                    gTemp /= r_a;
                    bTemp /= r_a;
                }
                
                if ( bd == 4 )
                {
                    wptr[0] = rTemp * 255.f;
                    wptr[1] = gTemp * 255.f;
                    wptr[2] = bTemp * 255.f;
                    wptr[3] = r_a * 255.f;
                }
                else
                {
                    wptr[0] = rTemp * ( 255.f * r_a );
                    wptr[1] = gTemp * ( 255.f * r_a );
                    wptr[2] = bTemp * ( 255.f * r_a );
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
                    //wptr[3] = 0xFF;
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
            return false;

        unsigned img_w  = src->w();
        unsigned img_h  = src->h();
        unsigned img_sz = img_w * img_h;

        uchar* ptr = (uchar*)src->data()[0];

        val = MIN( 1.0f, MAX( 0.0f, val ) );

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
        float col1g     = 0.0f;
        float col2g     = 0.0f;
        float randfv    = 0.0f;
        uchar fill_r    = 0;
        uchar fill_g    = 0;
        uchar fill_b    = 0;
        uchar fill_a    = 0;
        uchar ref_c1r   = ( col1 & 0xFF000000 ) >> 24;
        uchar ref_c1g   = ( col1 & 0x00FF0000 ) >> 16;
        uchar ref_c1b   = ( col1 & 0x0000FF00 ) >> 8;
        uchar ref_c1a   = alpha1;
        uchar ref_c2r   = ( col2 & 0xFF000000 ) >> 24;
        uchar ref_c2g   = ( col2 & 0x00FF0000 ) >> 16;
        uchar ref_c2b   = ( col2 & 0x0000FF00 ) >> 8;
        uchar ref_c2a   = alpha2;

        #pragma omp parellel for
        for ( int cy=0; cy<h; cy++ )
        {
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

        return new Fl_RGB_Image( buffer, w, h, d );
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
        float col1g     = 0.0f;
        float col2g     = 0.0f;
        float randfv    = 0.0f;
        uchar fill_r    = 0;
        uchar fill_g    = 0;
        uchar fill_b    = 0;
        uchar fill_a    = 0;
        uchar ref_c1r   = ( col1 & 0xFF000000 ) >> 24;
        uchar ref_c1g   = ( col1 & 0x00FF0000 ) >> 16;
        uchar ref_c1b   = ( col1 & 0x0000FF00 ) >> 8;
        uchar ref_c1a   = alpha1;
        uchar ref_c2r   = ( col2 & 0xFF000000 ) >> 24;
        uchar ref_c2g   = ( col2 & 0x00FF0000 ) >> 16;
        uchar ref_c2b   = ( col2 & 0x0000FF00 ) >> 8;
        uchar ref_c2a   = alpha2;

        #pragma omp parellel for private( cy )
        for ( int cx=0; cx<w; cx++ )
        {
            col1g = ( ( float( w*10 ) - float( cx*10 )  ) * downscale_f ) / 2559.0f;
            col2g = 1.0f - col1g;

            fill_r = ( (float)ref_c1r * col1g ) + ( (float)ref_c2r * col2g );
            fill_g = ( (float)ref_c1g * col1g ) + ( (float)ref_c2g * col2g );
            fill_b = ( (float)ref_c1b * col1g ) + ( (float)ref_c2b * col2g );
            if ( d > 3 )
            {
                fill_a = ( (float)ref_c1a * col1g ) + ( (float)ref_c2a * col2g );
            }

            for ( int cy=0; cy<h; cy++ )
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

        return new Fl_RGB_Image( buffer, w, h, d );
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

inline void fl_imgtk_dla_plot( Fl_RGB_Image* img, int x, int y, Fl_Color col, float br )
{
    float alpha = (float)( col & 0x000000FF ) / 255.f;
    uchar col_b = ( col & 0x0000FF00 ) >> 8;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_r = ( col & 0xFF000000 ) >> 24;

    int w = img->w();
    int h = img->h();
    int d = img->d();

    br *= alpha;

    if ( ( ( x >= 0 ) && ( y >= 0 ) )
        && ( ( x < w ) && ( y < h ) ) )
    {
        unsigned pos  = ( y * w + x ) * d;
        float revbr = 1.0f - br;

        uchar* ptrimg = (uchar*)img->data()[0];

        ptrimg[ pos + 0 ] = ( ptrimg[ pos + 0 ] * revbr ) + ( (float)col_r * br );
        ptrimg[ pos + 1 ] = ( ptrimg[ pos + 1 ] * revbr ) + ( (float)col_g * br );
        ptrimg[ pos + 2 ] = ( ptrimg[ pos + 2 ] * revbr ) + ( (float)col_b * br );
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

void fl_imgtk::draw_smooth_line( Fl_RGB_Image* img, int x1, int y1, int x2, int y2, Fl_Color col )
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

#define fl_imgtk_putpixel( _buff_,_x_,_w_,_y_,_d_,_r_,_g_,_b_,_a_ ) \
        uchar* _putptr_ = &_buff_[ ( ( _y_ * _w_ ) + _x_ ) * _d_ ];\
        if((uchar)_a_==0xFF)\
        {_putptr_[0]=_r_; _putptr_[1]=_g_; _putptr_[2]=_b_;}\
        else\
        {float _ar_=(float)_a_/255.f; float _rar_=1.f-_ar_;\
         _putptr_[0]=(uchar)((float)_putptr_[0]*_rar_) + ((float)_r_*_ar_);\
         _putptr_[1]=(uchar)((float)_putptr_[1]*_rar_) + ((float)_g_*_ar_);\
         _putptr_[2]=(uchar)((float)_putptr_[2]*_rar_) + ((float)_b_*_ar_); }
/*
inline void fl_imgtk_putpixel( uchar* _buff_, unsigned _x_, unsigned _w_, unsigned _y_, unsigned _d_, 
                               unsigned _r_, unsigned _g_, unsigned _b_, unsigned _a_ )
{
    uchar* _putptr_ = &_buff_[ ( ( _y_ * _w_ ) + _x_ ) * _d_ ];
    
    if((uchar)_a_==0xFF)
    {
        _putptr_[0]=_r_;
        _putptr_[1]=_g_;
        _putptr_[2]=_b_;
        
        if ( _d_ > 3 )
        {
            _putptr_[3] = (uchar)(( (float)_putptr_[3] / 255.f * (float)_a_/255.f ) * 255.f );
        }
    }
    else
    {
        float _ar_=(float)_a_ / 255.f; 
        float _rar_=1.f-_ar_;
        
        if ( _d_ < 4 )
        {
            _putptr_[0]=(uchar)((float)_putptr_[0]*_rar_) + ((float)_r_*_ar_);
            _putptr_[1]=(uchar)((float)_putptr_[1]*_rar_) + ((float)_g_*_ar_);
            _putptr_[2]=(uchar)((float)_putptr_[2]*_rar_) + ((float)_b_*_ar_);
        }
        else
        {
            _putptr_[0]=(uchar)((float)_putptr_[0]*_rar_) + ((float)_r_*_ar_);
            _putptr_[1]=(uchar)((float)_putptr_[1]*_rar_) + ((float)_g_*_ar_);
            _putptr_[2]=(uchar)((float)_putptr_[2]*_rar_) + ((float)_b_*_ar_);
            _putptr_[3]=(uchar)(((float)_putptr_[3]/255.f) * _ar_* 255.f);
        }
    }
}
*/
void fl_imgtk::draw_line( Fl_RGB_Image* img, int x1, int y1, int x2, int y2, Fl_Color col )
{
    if ( img == NULL )
    {
        return;
    }
    
    if ( img->d() < 3 )
    {
        return;
    }

    uchar* putbuff = (uchar*)img->data()[0];
    uchar col_b = ( col & 0x0000FF00 ) >> 8;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_r = ( col & 0xFF000000 ) >> 24;
    uchar col_a = ( col & 0x000000FF );
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

void fl_imgtk::draw_rect( Fl_RGB_Image* img, int x, int y, int w, int h, Fl_Color col )
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

void fl_imgtk::draw_fillrect( Fl_RGB_Image* img, int x, int y, int w, int h, Fl_Color col )
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
    uchar img_d = img->d();
    int   img_w = img->w();
    int   img_h = img->h();

    int cym = y+h;
    int cxm = x+w;

    if ( cxm > img_w )
        cxm = img_w - x;

    if ( cym > img_h )
        cym = img_h - y;

    #pragma omp parallel for
    for( int cy=y; cy<cym; cy++ )
    {
        for( int cx=x; cx<cxm; cx++ )
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

void fl_imgtk::draw_polygon( Fl_RGB_Image* img, const fl_imgtk::vecpoint* points, unsigned pointscnt, Fl_Color col )
{
    if ( img == NULL )
        return;
    
    if ( ( points == NULL ) || ( pointscnt < 3 ) )
        return;
    
    if ( img->d() < 3 )
        return;
    
    uchar col_r = ( col & 0x0000FF00 ) >> 8;
    uchar col_g = ( col & 0x00FF0000 ) >> 16;
    uchar col_b = ( col & 0xFF000000 ) >> 24;

    uchar* ptrimg = (uchar*)img->data()[0];
    
    const unsigned max_y = img->h();
    const unsigned max_x = img->w();

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

        unsigned node_x_sz = node_x.size();

        // sort nodes ..
        if ( node_x_sz > 1 )
        {
            sort( node_x.begin(),
                  node_x.begin() + node_x_sz,
                  fl_imgtk_sortcondition );

            for( unsigned ncnt=0; ncnt<node_x_sz; ncnt++ )
            {
                printf("%.2f, ", node_x[ncnt] );
            }
            printf("\n");

        }

        #pragma parallel for
        for( unsigned dcnt=0; dcnt<node_x_sz; dcnt+=2 )
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

                for( int cur_x=node_x[dcnt]; cur_x<=node_x[dcnt+1]; cur_x++ )
                {
                    int buffpos = ( max_x * cur_y + cur_x ) * img->d();
                    
                    ptrimg[ buffpos + 0 ] = col_r;
                    ptrimg[ buffpos + 1 ] = col_g;
                    ptrimg[ buffpos + 2 ] = col_b;
                }
            }
        }

        node_x.clear();

    } /// of for( y .. )

    img->uncache();
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
