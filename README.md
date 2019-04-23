# (lib)fl_imgtk

* A small library for FLTK image toolkit.
* Designed to use some useful effects in FLTK GUI.

## Lastest upate

* 2019-04-23-2
    1. New feature, draw_rect() and draw_fillrect() included.
	1. Version updated to 0.3.32.0

## Previous update

* 2019-04-23-1
    1. Little bug fixed in drawing simple line
	1. Version updated to 0.3.31.2

* 2019-04-22-1
    1. Little changes in drawing smooth line.
	1. removed non-affective OpenMP preprocessors.
	1. Version updated to 0.3.31.1

* 2018-03-26-1
    1. Fixed a bug of image merging.
    1. fl_drawonimage supports both of rgb and rgba images.
	1. Version updated to 0.3.31.0

* 2018-03-21-1
    1. Fixed a bug occurs black line when it turns image to 90 or 240.
	1. Version updated to 0.3.30.6

* 2018-02-02-1
    1. Fixed bugs in rotate free to use right degree in 0 to 359.
	1. Version updated to 0.3.30.5

* 2017-11-30-1
    1. Fixed bugs in draw_line() and draw_smoothline().
	1. Still left a bug in draw_smoothline() rounding error of y2 coordination.
	1. Version updated to 0.3.30.3

* 2017-11-14-1
    1. Fixed a bug in applyaplha_ex().
	1. VisualStudio 2015 project builds target for static build.
	1. Version updated to 0.3.30.2

* 2017-11-02-1
    1. Included a drawing function for single line (not aliased)
	1. Version updated to 0.3.30.1

* 2017-10-30-2
    1. Enhanced to OpenMP error for reinhard HDR, but running speed gone for slower.
	1. Version updated to 0.3.27.2
	
* 2017-10-30-1
    1. Fixed a bug of color order of makeanempty();
	1. Version updated to 0.3.27.1

* 2017-10-13-2
    1. Fixed SIGTRAP occurs when try remove Fl_RGB_Image buffer from got draw_currentwindow();
	1. It seems to a bug of FLTK, by fl_read_image().
	1. Version updated to 0.3.26.2

* 2017-10-13-1
    1. Changed sens bit depth in makegradation_h and v.
	1. Version updated to 0.3.26.1

* 2017-10-12-1
    1. Fixed a wrong method name, "brightbess_ex" to "brightness_ex".
	1. Added new function for "draw_currentwindow()".
	1. draw_smoothline() now controls with alpha channel.
	1. Version updated to 0.3.26.0

* 2017-10-11-1
    1. New method : createBMPmemory(), for create new BMP image from memory.
	1. Version updated to 0.3.25.1

* 2017-08-14-1
    1. New method : makeanempty(), for create new image.

* 2017-08-14-0
    1. Updated header for drawonimage method.
	1. drawonimage condition changed for alpha channel.

* 2017-08-10-0
    1. Fixed bugs in makegradation_h and v methods.
	1. Now makegradation_h/v() methods support alpha channel.
		
* 2017-07-27-0
    1. OpenMP applied drawing to 
	    - smooth line.
	    - polygon.

* 2017-07-26-0
    1. New functions included !
	    - drawing smooth line.
        - drawing polygon.		

* 2017-04-19-0
	1. Fixed wrong pointer addressing in
		- CLAHE();
		- CLAHE_ex();
	1. Added new effect belong to CLAHE.
	    - noire();
	    - noire_ex();

* 2017-04-18-0
	1. Added a feature
		- CLAHE();
		- CLAHE_ex();

* 2017-04-06-0
	1. Bug fixed
		- applyalpha_ex();

* 2017-04-05-0 
    1. New feature
        - invert();
        - invert_ex();

* 2017-03-31-0
    1. New features
        - makegradation_h();
        - makegradation_v();

* * 2017-03-30-0
    1. New features
        - FL_IMGTK, FL_IMGTK_VERSION definitions
        - fliphorizontal_ex();
        - flipvertical_ex();
        - gamma_ex();
        - brightbess_ex();
        - contrast_ex();
        - blurredimage_ex();
        - filtered_ex();
        - tonemapping_reinhard_ex();
        - tonemapping_drago_ex();
        - subtract();
        - sbutract_ex();
        - applyalpha_ex();

* * 2017-03-29-0
    1. Now supports many drawing image functions.
        - merge();
        - makealphamap();
        - applyalpha();
        - drawonimage();
        - 
* 2017-03-28-0
    1. Cropping image to a new Fl_RGB_Image.
    2. Fixed bugs on Tone mapping (HDR), now it runs well for Drago and Reinhard algorithms.

* 2017-03-24-0
    1. First commit.
    2. Supporting LLVM-gcc for Apple Mac.

## How to build ?

* Copy your right Makefile.xxxx to Makefile.
    - Ex) MinGW-W64 need do 
    ```cp Makefile.mingw Makefile```
* Then ```make```


## Currently supported

* OpenMP 3.0 or above with gcc, or supported compiler.
   - Some visual studio not process unsigned for() loop.
* Flip Fl_RGB_Image to horizental or vertical.
* Rotate Fl_RGB_Image to these :
    1. 90, 180, 270 degrees.
    2. Or free rotate
* Fast rescale with fl_smimg, with supporting these:
    1. Box(nearest)
    2. Bilinear
    3. Bicubic
    4. Lanczos
    5. B-Spline
* Draw Fl_Widget to Fl_RGB_Image.
* Draw Fl_Widget to blurred Fl_RGB_Image.
    - may useful to make background image.
* Tone mapping (HDR)
    1. DRAGO algorithm.
    2. REINHARD algorithm.
* Color CLAHE
* Kernel matrix filter ( 3x3 to NxN )
* Crop Fl_RGB_Image to a new Fl_RGB_Image.
* Merge two different Fl_RGB_Images to a new one.
* Make alpha map, and apply to Fl_RGB_Image.
* Draw Fl_RGB_Image to Fl_RGB_Image, it is same like Mix.
* Subtracts two Fl_RGB_Image.
* Make gradation image col1 to col2.

## Planned to next

* Channel mix
* White balance

## OpenMP supports

* fl_imgtk designed for using OpenMP 3.x or above,
* Some lower version of OpenMP may not work with unsigned loop.
* Apple Mac, X-Code may not supports OpenMP, recommend to use HPC-GCC to build library.
* To disable using OpenMP, remove -DUSING_OMP in makefile.

## Supported compiler:

* gcc (with OpenMP), or MinGW-W64
* llvm-gcc (Apple)
* (Maybe works on M$VC too, but no project file.)

## External licenses:

* FreeImage 3 (many codes belongs to this)
* FLTK License
