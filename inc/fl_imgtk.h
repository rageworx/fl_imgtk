#ifndef __FL_IMGTOOLKIT_H__
#define __FL_IMGTOOLKIT_H__

/*******************************************************************************
* fl_imgtk.H , version 2017-03-31-0
* =============================================================================
* A tool kit for basic FLTK image processing.
* (C) 2016-2017 Raphael Kim, Raph.K. ( rageworx or rage.kim @gmail.com )
* All rights reserved for MIT license.
*
* [ Disclaimer ]
* - Some codes belong to FreeImage 3 library, and modified to FLTK and fl_imgtk.
* - It follows FreeImage license and open as free.
*******************************************************************************/

#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_RGB_Image.H>

////////////////////////////////////////////////////////////////////////////////

#define FL_IMGTK

#define FL_IMGTK_VER_MJR    0
#define FL_IMGTK_VER_MNR    3
#define FL_IMGTK_VER_BLD    15
#define FL_IMGTK_VER_REV    0

#define FL_IMGTK_VERSION    ( FL_IMGTK_VER_MJR * 100000000 + \
                              FL_IMGTK_VER_MNR * 100000 + \
                              FL_IMGTK_VER_BLD * 1000 + \
                              FL_IMGTK_VER_REV )

////////////////////////////////////////////////////////////////////////////////

namespace fl_imgtk
{
    ////////////////////////////////////////////////////////////////////////////
    // enumerations, structures

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
    // Flip, Rotate

    Fl_RGB_Image* fliphorizontal( Fl_RGB_Image* img );
    bool          fliphorizontal_ex( Fl_RGB_Image* img );
    Fl_RGB_Image* flipvertical( Fl_RGB_Image* img );
    bool          flipvertical_ex( Fl_RGB_Image* img );

    Fl_RGB_Image* rotate90( Fl_RGB_Image* img );
    Fl_RGB_Image* rotate180( Fl_RGB_Image* img );
    Fl_RGB_Image* rotate270( Fl_RGB_Image* img );

    Fl_RGB_Image* rotatefree( Fl_RGB_Image* img, float deg );

    ////////////////////////////////////////////////////////////////////////////
    // Adjusting colors

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
    bool          gamma_ex( Fl_RGB_Image* img, double gamma );
    // perc => -100 ~ 100
    Fl_RGB_Image* brightness( Fl_RGB_Image* img, double perc );
    bool          brightbess_ex( Fl_RGB_Image* img, double perc );
    // perc => -100 ~ 100
    Fl_RGB_Image* contrast( Fl_RGB_Image* img, double perc );
    bool          contrast_ex(  Fl_RGB_Image* img, double perc );

    ////////////////////////////////////////////////////////////////////////////
    // Color conversion

    Fl_RGB_Image* invert( Fl_RGB_Image* img );
    bool          invert_ex( Fl_RGB_Image* img );

    ////////////////////////////////////////////////////////////////////////////
    // Rescale(resize) with fl_smimg

    Fl_RGB_Image* rescale( Fl_RGB_Image* img, unsigned w, unsigned h,
                           rescaletype rst = NONE );

    ////////////////////////////////////////////////////////////////////////////
    // Get/Drawing images

    Fl_RGB_Image* draw_widgetimage( Fl_Widget* w );
    Fl_RGB_Image* drawblurred_widgetimage( Fl_Widget* w, unsigned factor = 2 );
    Fl_RGB_Image* blurredimage( Fl_RGB_Image* src, unsigned factor = 2 );
    bool          blurredimage_ex( Fl_RGB_Image* src, unsigned factor = 2 );

    ////////////////////////////////////////////////////////////////////////////
    // Kernel matrix filter configurations
    /***
    * preset configs : blur, blurmore, sharpen, sharpenmore
    * null = empty one.
    ****/
    kfconfig* new_kfconfig( const char* preset );

    // Kernel matrix filter
    Fl_RGB_Image* filtered( Fl_RGB_Image* img, kfconfig* kfc );
    bool          filtered_ex( Fl_RGB_Image* img, kfconfig* kfc );

    ////////////////////////////////////////////////////////////////////////////
    // Tone mapping

	/***
	** tonemapping_reinhard
	** -------------------------------------------------------------------------
	** Apply the global/local tone mapping operator to a Fl_RGB_Image.
	** User parameters control intensity, contrast, and level of adaptation
	**
	** @param 'src': Input Fl_RGB_Image, depth requires at least 3.
	** @param 'intensity': Overall intensity in range [-8:8] : default to 0
	** @param 'contrast': Contrast in range [0.3:1) : default to 0
	** @param 'adaptation': Adaptation in range [0:1] : default to 1
	** @param 'color_correction': Color correction in range [0:1] : default to 0
	***/
    Fl_RGB_Image* tonemapping_reinhard( Fl_RGB_Image* src,
	                                    float intensity, float contrast,
										float adaptation = 1.0f,
										float color_correction = 0.0f );
    bool       tonemapping_reinhard_ex( Fl_RGB_Image* src,
	                                    float intensity, float contrast,
										float adaptation = 1.0f,
										float color_correction = 0.0f );

	/***
	** tonemapping_drago
	** --------------------------------------------------------------------------
	** Apply the Adaptive Logarithmic Mapping operator to a HDR image.
	**
	** @param 'src': Input Fl_RGB_Image, depth requires at least 3.
	** @param 'gamma': Gamma correction (gamma > 0). 1 means no correction,
	**                                               2.2 in the original paper.
	** @param 'exposure': Exposure parameter
	**                        (0 means no correction, 0 in the original paper)
	***/
    Fl_RGB_Image* tonemapping_drago( Fl_RGB_Image* src,
                                     float gamma = 1.0f ,
                                     float exposure = 0.0f );
    bool       tonemapping_drago_ex( Fl_RGB_Image* src,
                                     float gamma = 1.0f ,
                                     float exposure = 0.0f );


    ////////////////////////////////////////////////////////////////////////////
    // More functions ...

	typedef struct
	{
	    int src1putx;
	    int src1puty;
	    int src2putx;
	    int src2puty;

	    float src1ratio;    /// 0.0 to 1.0
	    float src2ratio;    /// 0.0 to 1.0
	                        /// src1ratio + src2ratio must be 1.0
        bool  autoexpand;   /// expands maximum image size to bigger image.
	}mergeconfig;

	// Crop image to a new Fl_RGB_Image.
    Fl_RGB_Image* crop( Fl_RGB_Image* src,
                        unsigned sx, unsigned sy, unsigned w, unsigned h );

    // Merge two different image to a new Image.
	Fl_RGB_Image* merge( Fl_RGB_Image* src1, Fl_RGB_Image* src2,
                         mergeconfig* cfg = NULL );

    // Subtract twro different image src1 - ( src2 * sr ) = result.
    Fl_RGB_Image* subtract( Fl_RGB_Image* src1, Fl_RGB_Image* src2,
                            int px, int py, float sr = 1.0f );
    // Subtract src1 - ( src2 * sr ) = src1.
    bool          sbutract_ex( Fl_RGB_Image* src1, Fl_RGB_Image* src2,
                               int px, int py, float sr = 1.0f );

    // Makes alpha channel map.
    unsigned makealphamap( uchar* &amap, Fl_RGB_Image* src, float val = 1.0f );
    unsigned makealphamap( uchar* &amap, unsigned w, unsigned h, uchar val  = 255 );

    // Apply alpha channel map to Fl_RGB_Image.
    // return Fl_RGB_Image always became to 4 depth image.
    Fl_RGB_Image* applyalpha( Fl_RGB_Image* src, float val = 1.0f );
    Fl_RGB_Image* applyalpha( Fl_RGB_Image* src,
                              uchar* alphamap = NULL , unsigned amsz = 0);
    bool          applyalpha_ex( Fl_RGB_Image* src, float val = 1.0f );

    // Draws image to image with alpha value.
    bool drawonimage( Fl_RGB_Image* bgimg, Fl_RGB_Image* img,
                      int x, int y, float alpha = 1.0f );

    // Make a RGB image from gradation col1 to col2.
    // Use col1, col2 to 0xRRGGBB00 channel alignment, not Fl_Color !
    Fl_RGB_Image* makegradation_h( unsigned w, unsigned h,
                                   ulong col1, ulong col2, bool dither );
    Fl_RGB_Image* makegradation_v( unsigned w, unsigned h,
                                   ulong col1, ulong col2, bool dither );

    ////////////////////////////////////////////////////////////////////////////

    void discard_user_rgb_image( Fl_RGB_Image* &img );
    void discard_kfconfig( kfconfig* &kfc ); /// Not Kentucky Fried Chicken !
};

#endif /// of __FL_IMGTOOLKIT_H__
