Requires:

QGLViewer (known to work with version 2.2.1)
http://artis.inrialpes.fr/Software/QGLViewer/

lib3ds
http://lib3ds.sourceforge.net/
(
You can use the 1.2.0 release if working on fedora core 3.
However fedora core 4 uses gcc4 and the 1.2.0 release of lib3ds
contains a bug which causes model loading to fail silently if compiled 
with gcc4. This bug is fixed in the 1.3.0 release which is
available by grabbing the lib3ds sourcecode via cvs and compiling
from that. Use the commands from README-dist to get a configure script.
)
