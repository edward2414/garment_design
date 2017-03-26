//
//  Filename         : offscreen.h
//  Author           : Cyril Soler
//  Purpose          : An Offscreen renderer.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  OFFSCREEN_H
# define OFFSCREEN_H

# include <stdlib.h>
# include <stdio.h>
# include <iostream>
# include <malloc.h>
# include <cmath>
# include <limits.h>
# include <sys/types.h>
# include <sys/times.h>

# include <X11/Xlib.h>
# include <GL/gl.h>
# include <GL/glx.h>
# include <GL/glu.h>

// I put this here to avoid conflicts between Qt and
// X11 in the definition of INT32

namespace offscreen {

  /*! \namespace offscreen
   *  \bugs attention a la restauration de context quand il n'avait rien avant
   *  l'OffScreenArea !!
   */

  static bool OutOfMemory = false ;
  static XErrorHandler oldHandler = NULL ;

  static int myXErrorHandler (Display *_d, XErrorEvent *_xee)
  {
    OutOfMemory = True;
    if(oldHandler)
      {
	return (*oldHandler)(_d,_xee);
      }
    else
      {
	return false;
      }

    return 0;
  }

  class OffScreenArea 
  {
  public:
    static const int UNKNOWN_OFFSCREEN_TYPE = 0 ;
    static const int PIXMAP_OFFSCREEN_TYPE = 1 ;
    static const int PBUFFER_OFFSCREEN_TYPE = 2 ;
    static const int KEEP_PIXMAP = 1 ;

    OffScreenArea  (int type = UNKNOWN_OFFSCREEN_TYPE, GLXContext shareCtx = NULL)
    {
      DefaultType = type ;
      i_Screen = -1;
      i_Height = 0;
      i_Width  = 0;
			
      i_OffScreenAreaType = UNKNOWN_OFFSCREEN_TYPE ;
      i_GLContext = NULL;
      i_shareContext = shareCtx;
      i_Drawable  = GLXDrawable(0);
      i_XPix      = Pixmap(0);
      i_pXVisual  = NULL;
      i_pDisplay  = NULL;

      SetDisplayAndScreen((Display *)NULL,(int)-1) ;
    }

    ~OffScreenArea ()
    {
      DestroyOffScreenArea();
    }

    GLXContext GetGLXContext() {
      return i_GLContext ;
    }
	
    /*
     * 0 : cannot allocate
     * 1 : allocation done
     * FIRST : try to allocate PixelBuffer (Single, double buffer)
     * SECOND : try to allocate GLXPixmap if PixelBuffer is not available
     */
    /*
     * Here, we try to allocate an OffScreen Area
     * first, with PBuffer (if this GLX_ext. is
     * available and we can create one)
     * second, with a pixmap 
     */
    int AllocateOffScreenArea(int width,int height)
    {
      // sanity check. Did assigning a display work?
      if(i_pDisplay == NULL)
      {
        std::cerr << "X Problem: Could not connect to a display" << std::endl;
        return 0;
      }
      
      
      save_GLContext = glXGetCurrentContext();
      save_pDisplay  = glXGetCurrentDisplay();
      save_Drawable  = glXGetCurrentDrawable();

      int AlreadyThere = 0;
# ifdef A_VIRER
      static int AlreadyThere = 0;

      if ( ( width != i_Width ) || ( height != i_Height ) )
	{
	  AlreadyThere = 0;
	  DestroyOffScreenArea();
	}
# endif
      if(!AlreadyThere)
	{
	  AlreadyThere = 1;

	  /** Before to use a Pixmap, we try with pbuffer **/	
	  if(TryPBuffer(false,width,height))
	    {
# ifdef DEBUG
	      fprintf(stderr, "Using single-buffer PBuffer for off-screen rendering.\n") ;
# endif
	      return true ;
	    }

	  fprintf(stderr, "Cannot use a single-buffer PBuffer, trying double buffer.\n") ;

	  if(TryPBuffer(true,width,height))
	    {
# ifdef DEBUG
	      fprintf(stderr, "Using double-buffer PBuffer for off-screen rendering.\n") ;
# endif
	      return true ;
	    }
# ifdef DEBUG
	  fprintf(stderr, "Warning : cannot create a PBuffer, trying Pixmap\n");
# endif
	  if(TryPixmap(width,height))
	    {
# ifdef DEBUG
	      fprintf(stderr, "Notice  : using Pixmap for offScreen rendering\n");
# endif
	      return true ;
	    }
# ifdef DEBUG
	  fprintf (stderr, "Warning : cannot create a Pixmap\n");
# endif

	  return false;
	}
      return true ;
    }

    void MakeCurrent() 
    {
      glXMakeCurrent(i_pDisplay, i_Drawable, i_GLContext);
    }
	
  protected:

    int DefaultType ;
    int i_OffScreenAreaType;
    int i_Height;
    int i_Width;
	
    GLXContext    save_GLContext;
    Display     * save_pDisplay;
    GLXDrawable   save_Drawable;

    Display     * i_pDisplay;
    int           i_Screen;
    GLXContext    i_GLContext;
    GLXContext    i_shareContext;
    GLXDrawable   i_Drawable;
    Pixmap        i_XPix;
    XVisualInfo * i_pXVisual;
	
    /* 
     * Define Display and screen
     * IF display == NULL THEN try to open default Display
     * IF screenNumber is < 0 THEN take default Screen
     */
    void SetDisplayAndScreen ( Display *pDisplay , int Screen )
    {
      //JDW The code was designed with this call in the class constructor. 
      //JDW. But as these calls can fail this code shouldn't be in the constructor.
      //JDW These warning were added to help debugging of problems in this case. 
      //JDW usually due to no DISPLAY variable set
      
      if ( pDisplay == NULL )
      {
	     i_pDisplay = XOpenDisplay ( NULL );
      }
      else
	   i_pDisplay = pDisplay;
      
      if(i_pDisplay == NULL)
      {
        std::cerr << "X Problem, could not connect to a display. Check DISPLAY env. var." << std::endl;
      }
			
      if ( Screen < 0 )
        i_Screen = DefaultScreen ( i_pDisplay );
      else
        i_Screen = Screen;
      
      if(i_Screen<0)
      {
        std::cerr << "X Problem, could not assign a screen. Got screen number: " << i_Screen << std::endl;
      }
    }
		
    // 
    // Creates a PBuffer
    // 
    // <Return Values>
    // 0 : failure
    // 1 : succeed
    //  

    bool CreatePBuffer (unsigned int width, unsigned int height , int * pAttribList)
    {
# ifdef DEBUG	
      int error = 0 ;
      while((error = glGetError()) > 0)
	std::cerr << "GLError " << (void *)error << " encountered." << std::endl ;
# endif
      GLXFBConfig *pfbConfigs;
      int nbConfigs;
      static int pbAttribs[] = { GLX_LARGEST_PBUFFER, true,
				 GLX_PRESERVED_CONTENTS, true,
				 GLX_PBUFFER_WIDTH,0,
				 GLX_PBUFFER_HEIGHT,0,
				 None };
			
      pbAttribs[5] = width ;
      pbAttribs[7] = height ;

      // Looks for a config that matches pAttribList
      pfbConfigs = glXChooseFBConfig(i_pDisplay,i_Screen,pAttribList,&nbConfigs) ;
# ifdef DEBUG
      std::cout << nbConfigs << " found for pbuffer." << std::endl ;
# endif

      if(pfbConfigs == NULL)
      {
        std::cerr << "Could not get FBConfig" << std::endl;
          return false ;
      }

      i_pXVisual = glXGetVisualFromFBConfig(i_pDisplay,pfbConfigs[0]);
      if(i_pXVisual == NULL)
      {
        std::cerr << "Could not get required type of display" << std::endl;
        return false ;
      }
      
      
      
      i_OffScreenAreaType = PBUFFER_OFFSCREEN_TYPE;

      // Sets current error handler
      OutOfMemory = False;
      oldHandler = XSetErrorHandler( myXErrorHandler );

      i_Drawable = glXCreatePbuffer(i_pDisplay,pfbConfigs[0],pbAttribs);
					
      if(i_Drawable == 0)
	{
	  i_pXVisual = NULL;
	  return false ;
	}
      unsigned int w=0,h=0;
      glXQueryDrawable(i_pDisplay,i_Drawable,GLX_WIDTH,&w) ;
      glXQueryDrawable(i_pDisplay,i_Drawable,GLX_HEIGHT,&h) ;

      if((w != width)||(h != height))
	{
# ifdef DEBUG
	  std::cerr << "Could not allocate Pbuffer. Only size " << w << "x" << h << " found." << std::endl ;
# endif
	  return false ;
	}
# ifdef DEBUG
      else
	std::cerr << "Could allocate Pbuffer. Size " << w << "x" << h << " found." << std::endl ;
# endif
# ifdef DEBUG	
      while((error = glGetError()) > 0)
	std::cerr << "GLError " << (void *)error << " encountered." << std::endl ;
# endif
      // now create GLXContext

      if((i_GLContext = glXCreateContext(i_pDisplay,i_pXVisual,i_shareContext,true)) == NULL)
	{
	  DestroyOffScreenArea() ;
	  return false ;
	}

      /* Restore original X error handler */
      (void) XSetErrorHandler( oldHandler );
		
      if(!OutOfMemory)
	{	
	  i_Height = height;
	  i_Width  = width;
			
	  return true ;
	}
      else
	return false ;
    }

    // 
    // Creates a Pixmap
    // 
    // <Return Values>
    // false : failure
    // true  : succeed
    // 

    bool CreatePixmap (int width, int height , int * pAttribList)
    {
      int depth;
      int totdepth=0;
      XErrorHandler oldHandler;
      XVisualInfo * pvisP;
			
      pvisP = glXChooseVisual ( i_pDisplay, i_Screen , pAttribList);

      if ( pvisP == NULL)
	{
	  fprintf( stderr , "Warning : no 24-bit true color visual available\n" );
	  return false ;
	}
			
      OutOfMemory = False;
      oldHandler = XSetErrorHandler(myXErrorHandler);
      if(i_XPix == Pixmap(NULL))
	{
	  depth = 0;
	  for (unsigned int i=0,j=0; (pAttribList[i] != None) && (j<3) ; i++ )
	    {
	      switch ( pAttribList[i] )
		{
		case GLX_RED_SIZE: 	glXGetConfig(i_pDisplay,pvisP,GLX_RED_SIZE,&depth) ;
		  totdepth += depth ;
		  i++ ;
		  j++ ;
		  break;
		
		case GLX_GREEN_SIZE: glXGetConfig(i_pDisplay,pvisP,GLX_GREEN_SIZE,&depth) ;
		  totdepth += depth ;
		  i++ ;
		  j++ ;
		  break;
		
		case GLX_BLUE_SIZE: 	glXGetConfig(i_pDisplay,pvisP,GLX_BLUE_SIZE,&depth) ;
		  totdepth += depth ;
		  i++ ;
		  j++ ;
		  break;

		case GLX_ALPHA_SIZE: 	glXGetConfig(i_pDisplay,pvisP,GLX_ALPHA_SIZE,&depth) ;
		  totdepth += depth ;
		  i++ ;
		  j++ ;
		  break;

		default:
		  break;
		}
	    }
			    
	  fprintf(stderr,"%d bits color buffer found\n",depth) ;
	  i_XPix = XCreatePixmap(i_pDisplay,RootWindow (i_pDisplay,0),width,height,totdepth);
	  XSync(i_pDisplay,False);
	  if(OutOfMemory)
	    {
	      i_XPix = Pixmap(0);
	      XSetErrorHandler(oldHandler);
	      oldHandler = NULL ;
	      fprintf(stderr,"Warning : could not allocate Pixmap\n");
	      return false ;
	    }
	}
		
      // Perhaps should we verify th type of Area (Pixmap) ?
      if ( i_Drawable == GLXDrawable(NULL) )
	{
	  // i_Drawable = i_XPix;
	  i_Drawable = glXCreateGLXPixmap ( i_pDisplay , pvisP , i_XPix );
	  XSync ( i_pDisplay , False );
	  if(OutOfMemory)
	    { 
	      i_Drawable = GLXDrawable(0); 
	      DestroyOffScreenArea();
	      fprintf ( stderr , "Warning : could not allocate GLX Pixmap\n"); 
	      return false ;
	    } 
	  else
	    {
	      if(i_GLContext != NULL) 
		{
		  glXDestroyContext ( i_pDisplay , i_GLContext );
		  i_GLContext = NULL;
		}
	      if((i_GLContext = glXCreateContext(i_pDisplay,pvisP,NULL,GL_FALSE)) == NULL)
		{
		  DestroyOffScreenArea();
		  fprintf(stderr, "Warning : could not create rendering context");
		}
	    }
	}
      XSetErrorHandler(oldHandler);
			
      i_pXVisual = (i_Drawable != GLXDrawable(NULL) ? pvisP : NULL);
			
      i_Height = height;
      i_Width  = width;
			
      if(i_Drawable != GLXDrawable(NULL))
	{
	  i_OffScreenAreaType = PIXMAP_OFFSCREEN_TYPE;
	  return true ;
	}
			
      return false ;
    }

    bool TryPixmap(int width,int height)
    {
      int attrList[30];
      int n = 0;
		
      attrList[n++] = GLX_RED_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_GREEN_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_BLUE_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_ALPHA_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_RGBA;
      attrList[n++] = GLX_DEPTH_SIZE;
      attrList[n++] = 16;
      attrList[n++] = GLX_STENCIL_SIZE;
      attrList[n++] = 1;
      attrList[n++] = None;
		
      return CreatePixmap(width,height,attrList) ;
    }

    bool TryPBuffer(bool double_buffer,int width,int height)
    {
      int attrList[30];
      int n = 0;
		
      attrList[n++] = GLX_RENDER_TYPE;
      attrList[n++] = GLX_RGBA_BIT;
      attrList[n++] = GLX_DRAWABLE_TYPE;
      attrList[n++] = GLX_PBUFFER_BIT;
      attrList[n++] = GLX_RED_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_GREEN_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_BLUE_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_ALPHA_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_DEPTH_SIZE;
      attrList[n++] = 8;
      attrList[n++] = GLX_DOUBLEBUFFER;
      attrList[n++] = double_buffer;
      attrList[n++] = None;
		
      return CreatePBuffer(width,height,attrList) ;
    }

    void DestroyOffScreenArea()
    {
      glXMakeCurrent(save_pDisplay,save_Drawable,save_GLContext);

      switch ( i_OffScreenAreaType )
	{
	case PIXMAP_OFFSCREEN_TYPE : 	if(i_Drawable != 0)
	  glXDestroyGLXPixmap(i_pDisplay,i_Drawable);
					
	  if (i_XPix != 0)
	    {
	      XFreePixmap (i_pDisplay,i_XPix);
	      i_XPix = 0;
	    }
	  break;
	case PBUFFER_OFFSCREEN_TYPE :	if(i_Drawable != 0)
	  glXDestroyPbuffer(i_pDisplay,i_Drawable);
	  break;
	default:
	  ;
	}
		
      if (i_GLContext != NULL)
	{
	  glXDestroyContext(i_pDisplay,i_GLContext);
	  i_GLContext = NULL;
	}
		  
      i_Drawable = 0;
      i_OffScreenAreaType = UNKNOWN_OFFSCREEN_TYPE;
    }

  };

} // end of namespace offscreen

#endif // OFFSCREEN_H
