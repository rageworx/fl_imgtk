# fl_imgtk

* * A small library for FLTK image toolkit
* Designed to use some useful effects in FLTK GUI.

## Lastest upate

* First commit.


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
* Tone mapping (HDR) - Experimental
    1. DRAGO algorithm.
    2. REINHARD algorithm.
* Kernel matrix filter ( 3x3 to NxN )

## Planned to next

* Crop image, or copy image to image, or merge images
* Crop image from polygonal vectors
* Apply or remove alpha channel
* Mix images

## Supported compiler:

* gcc (with OpenMP), or MinGW-W64


## External license:

* FreeImage 3 (many codes belongs to this)
* FLTK License