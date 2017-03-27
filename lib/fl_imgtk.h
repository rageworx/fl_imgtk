#ifndef __FL_IMGTOOLKIT_H__
#define __FL_IMGTOOLKIT_H__

/*******************************************************************************
* fl_imgtk.H , version 2017-03-24-0
* =============================================================================
* A tool kit for basic FLTK image processing.
* (C) 2016-2017 Raphael Kim, Raph.K. ( rageworx or rage.kim @gmail.com )
* All rights reserved for MIT license.
*******************************************************************************/

#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_RGB_Image.H>

namespace fl_imgtk
{
    typedef enum
    {
        NONE = 0,
        BILINEAR,
        BICUBIC,
        LANCZOS,
        BSPLINE
    }rescaletype;

    typedef struct
    {
        uchar   w;      /// Width
        uchar   h;      /// Height
        float   f;      /// Factor
        float   b;      /// Bias
        uchar   msz;    /// Matrix size
        float*  m;      /// Matrix array (dynamical allocation)
    }kfconfig; /// Kernel Filter Configuration.

    ////////////////////////////////////////////////////////////////////////////

    Fl_RGB_Image* fliphorizontal( Fl_RGB_Image* img );
    Fl_RGB_Image* flipvertical( Fl_RGB_Image* img );

    Fl_RGB_Image* rotate90( Fl_RGB_Image* img );
    Fl_RGB_Image* rotate180( Fl_RGB_Image* img );
    Fl_RGB_Image* rotate270( Fl_RGB_Image* img );

    Fl_RGB_Image* rotatefree( Fl_RGB_Image* img, float deg );

    /***
    * these following methods: Gamma, Brightness, Contrast
    * are motivated from FreeImageToolkit open source.
    * --  http://freeimage.sourceforge.net/
    *                      ,by Following FreeImage License.
    *
    * reprogrammed by Raph.K. (rageworx@gmail.com).
    ***/
    // Gamma default = 1.0
    Fl_RGB_Image* gamma( Fl_RGB_Image* img, double gamma );
    // perc => -100 ~ 100
    Fl_RGB_Image* brightness( Fl_RGB_Image* img, double perc );
    // perc => -100 ~ 100
    Fl_RGB_Image* contrast( Fl_RGB_Image* img, double perc );

    // Kernel matrix filter
    Fl_RGB_Image* filtered( Fl_RGB_Image* img, kfconfig* kfc );

    // Rescale(resize) with fl_smimg
    Fl_RGB_Image* rescale( Fl_RGB_Image* img, unsigned w, unsigned h, rescaletype rst = NONE );

    ////////////////////////////////////////////////////////////////////////////

    Fl_RGB_Image* draw_widgetimage( Fl_Widget* w );
    Fl_RGB_Image* drawblurred_widgetimage( Fl_Widget* w, unsigned factor = 2);

    /***
    * preset configs : blur, blurmore, sharpen, sharpenmore
    * null = empty one.
    ****/
    kfconfig* new_kfconfig( const char* preset );

    ////////////////////////////////////////////////////////////////////////////
    // Tone mapping

    Fl_RGB_Image* tonemapping_reinhard( Fl_RGB_Image* src, float intensity, float contrast, float adaptation = 1.0f, float color_correction = 0.0f );
    Fl_RGB_Image* tonemapping_drago( Fl_RGB_Image* src, float gamma = 1.0f , float exposure = 0.0f );
    ////////////////////////////////////////////////////////////////////////////

    void discard_user_rgb_image( Fl_RGB_Image* &img );
    void discard_kfconfig( kfconfig* &kfc ); /// Not Kentucky Fried Chicken !
};

#endif /// of __FL_IMGTOOLKIT_H__
