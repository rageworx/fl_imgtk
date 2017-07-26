# (lib)fl_imgtk

* A small library for FLTK image toolkit.
* Designed to use some useful effects in FLTK GUI.

## Lastest upate

* 2017-07-26-0
    1. New functions included !
	    - drawing smooth line.
        - drawing polygon.		

## Previous update

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
