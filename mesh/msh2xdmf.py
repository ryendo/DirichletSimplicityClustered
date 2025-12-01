#!/Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12

import meshio
msh = meshio.read("/tmp/temp.msh")
meshio.write("/tmp/temp.xdmf", msh) # Writes to XDMF, an XML-based format