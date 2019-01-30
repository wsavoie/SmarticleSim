###########################################################################
## Makefile generated for MATLAB file/project 'generatePackingFromSimDatPLANE'. 
## 
## Makefile     : generatePackingFromSimDatPLANE_rtw.mk
## Generated on : Tue Dec 18 12:31:23 2018
## MATLAB Coder version: 4.0 (R2018a)
## 
## Build Info:
## 
## Final product: $(RELATIVE_PATH_TO_ANCHOR)\generatePackingFromSimDatPLANE.lib
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPUTER                Computer type. See the MATLAB "computer" command.
# CMD_FILE                Command file

PRODUCT_NAME              = generatePackingFromSimDatPLANE
MAKEFILE                  = generatePackingFromSimDatPLANE_rtw.mk
COMPUTER                  = PCWIN64
MATLAB_ROOT               = C:\PROGRA~1\MATLAB\R2018a
MATLAB_BIN                = C:\PROGRA~1\MATLAB\R2018a\bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)\win64
MASTER_ANCHOR_DIR         = 
START_DIR                 = D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\newSimScripts
ARCH                      = win64
RELATIVE_PATH_TO_ANCHOR   = .
CMD_FILE                  = generatePackingFromSimDatPLANE_rtw.rsp
C_STANDARD_OPTS           = 
CPP_STANDARD_OPTS         = 
NODEBUG                   = 1

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          Microsoft Visual C++ 2017 v15.0 | nmake (64-bit Windows)
# Supported Version(s):    15.0
# ToolchainInfo Version:   R2018a
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# NODEBUG
# cvarsdll
# cvarsmt
# conlibsmt
# ldebug
# conflags
# cflags

#-----------
# MACROS
#-----------

MEX_OPTS_FILE       = $(MATLAB_ROOT)\bin\$(ARCH)\mexopts\msvc2017.xml
MW_EXTERNLIB_DIR    = $(MATLAB_ROOT)\extern\lib\win64\microsoft
MW_LIB_DIR          = $(MATLAB_ROOT)\lib\win64
MEX_ARCH            = -win64
CPU                 = AMD64
APPVER              = 5.02
CVARSFLAG           = $(cvarsmt)
CFLAGS_ADDITIONAL   = -D_CRT_SECURE_NO_WARNINGS
CPPFLAGS_ADDITIONAL = -EHs -D_CRT_SECURE_NO_WARNINGS
LIBS_TOOLCHAIN      = $(conlibs)

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: Microsoft Visual C Compiler
CC = cl

# Linker: Microsoft Visual C Linker
LD = link

# C++ Compiler: Microsoft Visual C++ Compiler
CPP = cl

# C++ Linker: Microsoft Visual C++ Linker
CPP_LD = link

# Archiver: Microsoft Visual C/C++ Archiver
AR = lib

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_ARCH_BIN)
MEX = "$(MEX_PATH)\mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: NMAKE Utility
MAKE = nmake


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -Zi
C_OUTPUT_FLAG       = -Fo
LDDEBUG             = /DEBUG
OUTPUT_FLAG         = -out:
CPPDEBUG            = -Zi
CPP_OUTPUT_FLAG     = -Fo
CPPLDDEBUG          = /DEBUG
OUTPUT_FLAG         = -out:
ARDEBUG             =
STATICLIB_OUTPUT_FLAG = -out:
MEX_DEBUG           = -g
RM                  = @del
ECHO                = @echo
MV                  = @ren
RUN                 = @cmd /C

#----------------------------------------
# "Faster Builds" Build Configuration
#----------------------------------------

ARFLAGS              = /nologo
CFLAGS               = $(cflags) $(CVARSFLAG) $(CFLAGS_ADDITIONAL) \
                       /Od /Oy-
CPPFLAGS             = /TP $(cflags) $(CVARSFLAG) $(CPPFLAGS_ADDITIONAL) \
                       /Od /Oy-
CPP_LDFLAGS          = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN)
CPP_SHAREDLIB_LDFLAGS  = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN) \
                         -dll -def:$(DEF_FILE)
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN)
MEX_CPPFLAGS         = -R2018a $(MEX_ARCH) OPTIMFLAGS="/Od /Oy- $(MDFLAG) $(DEFINES)" $(MEX_OPTS_FLAG)
MEX_CPPLDFLAGS       =
MEX_CFLAGS           = -R2018a $(MEX_ARCH) OPTIMFLAGS="/Od /Oy- $(MDFLAG) $(DEFINES)" $(MEX_OPTS_FLAG)
MEX_LDFLAGS          = LDFLAGS=='$$LDFLAGS'
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN) \
                       -dll -def:$(DEF_FILE)

#--------------------
# File extensions
#--------------------

H_EXT               = .h
OBJ_EXT             = .obj
C_EXT               = .c
EXE_EXT             = .exe
SHAREDLIB_EXT       = .dll
HPP_EXT             = .hpp
OBJ_EXT             = .obj
CPP_EXT             = .cpp
EXE_EXT             = .exe
SHAREDLIB_EXT       = .dll
STATICLIB_EXT       = .lib
MEX_EXT             = .mexw64
MAKE_EXT            = .mk


###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = $(RELATIVE_PATH_TO_ANCHOR)\generatePackingFromSimDatPLANE.lib
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE;$(START_DIR);$(MATLAB_ROOT)\extern\include;$(MATLAB_ROOT)\simulink\include;$(MATLAB_ROOT)\rtw\c\src;$(MATLAB_ROOT)\rtw\c\src\ext_mode\common;$(MATLAB_ROOT)\rtw\c\ert

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_STANDARD = -DMODEL=generatePackingFromSimDatPLANE -DHAVESTDIO -DUSE_RTMODEL

DEFINES = $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\generatePackingFromSimDatPLANE_data.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\generatePackingFromSimDatPLANE_initialize.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\generatePackingFromSimDatPLANE_terminate.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\generatePackingFromSimDatPLANE.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\quat2rotm.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\sum.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\bsxfun.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\pdist.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\intersectRayPolygon3d.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\createPlane.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\normalizeVector3d.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\crossProduct3d.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\projPointOnPlane.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\vectorNorm3d.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\dot.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\splitPolygons.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\polygonArea.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\inpolygon.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\distancePointEdge3d.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\generatePackingFromSimDatPLANE_emxutil.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\generatePackingFromSimDatPLANE_emxAPI.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\rt_nonfinite.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\rtGetNaN.cpp $(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE\rtGetInf.cpp

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = generatePackingFromSimDatPLANE_data.obj generatePackingFromSimDatPLANE_initialize.obj generatePackingFromSimDatPLANE_terminate.obj generatePackingFromSimDatPLANE.obj quat2rotm.obj sum.obj bsxfun.obj pdist.obj intersectRayPolygon3d.obj createPlane.obj normalizeVector3d.obj crossProduct3d.obj projPointOnPlane.obj vectorNorm3d.obj dot.obj splitPolygons.obj polygonArea.obj inpolygon.obj distancePointEdge3d.obj generatePackingFromSimDatPLANE_emxutil.obj generatePackingFromSimDatPLANE_emxAPI.obj rt_nonfinite.obj rtGetNaN.obj rtGetInf.obj

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = 

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = 

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_ = /openmp
CFLAGS_BASIC = $(DEFINES) 

CFLAGS = $(CFLAGS) $(CFLAGS_) $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_ = /openmp
CPPFLAGS_BASIC = $(DEFINES) 

CPPFLAGS = $(CPPFLAGS) $(CPPFLAGS_) $(CPPFLAGS_BASIC)

#---------------
# C++ Linker
#---------------

CPP_LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

CPP_LDFLAGS = $(CPP_LDFLAGS) $(CPP_LDFLAGS_)

#------------------------------
# C++ Shared Library Linker
#------------------------------

CPP_SHAREDLIB_LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS) $(CPP_SHAREDLIB_LDFLAGS_)

#-----------
# Linker
#-----------

LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

LDFLAGS = $(LDFLAGS) $(LDFLAGS_)

#--------------------------
# Shared Library Linker
#--------------------------

SHAREDLIB_LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS) $(SHAREDLIB_LDFLAGS_)

###########################################################################
## INLINED COMMANDS
###########################################################################


!include $(MATLAB_ROOT)\rtw\c\tools\vcdefs.mak


###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute set_environment_variables


all : build
	@cmd /C "@echo ### Successfully generated all binary outputs."


build : set_environment_variables prebuild $(PRODUCT)


prebuild : 


download : build


execute : download


set_environment_variables : 
	@set INCLUDE=$(INCLUDES);$(INCLUDE)
	@set LIB=$(LIB)


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	@cmd /C "@echo ### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS) -out:$(PRODUCT) @$(CMD_FILE)
	@cmd /C "@echo ### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\generatePackingFromSimDatPLANE}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\rtw\c\src}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\rtw\c\src}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : $(MAKEFILE) rtw_proj.tmw


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@cmd /C "@echo ### PRODUCT = $(PRODUCT)"
	@cmd /C "@echo ### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@cmd /C "@echo ### BUILD_TYPE = $(BUILD_TYPE)"
	@cmd /C "@echo ### INCLUDES = $(INCLUDES)"
	@cmd /C "@echo ### DEFINES = $(DEFINES)"
	@cmd /C "@echo ### ALL_SRCS = $(ALL_SRCS)"
	@cmd /C "@echo ### ALL_OBJS = $(ALL_OBJS)"
	@cmd /C "@echo ### LIBS = $(LIBS)"
	@cmd /C "@echo ### MODELREF_LIBS = $(MODELREF_LIBS)"
	@cmd /C "@echo ### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@cmd /C "@echo ### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@cmd /C "@echo ### CFLAGS = $(CFLAGS)"
	@cmd /C "@echo ### LDFLAGS = $(LDFLAGS)"
	@cmd /C "@echo ### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@cmd /C "@echo ### CPPFLAGS = $(CPPFLAGS)"
	@cmd /C "@echo ### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@cmd /C "@echo ### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@cmd /C "@echo ### ARFLAGS = $(ARFLAGS)"
	@cmd /C "@echo ### MEX_CFLAGS = $(MEX_CFLAGS)"
	@cmd /C "@echo ### MEX_CPPFLAGS = $(MEX_CPPFLAGS)"
	@cmd /C "@echo ### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@cmd /C "@echo ### MEX_CPPLDFLAGS = $(MEX_CPPLDFLAGS)"
	@cmd /C "@echo ### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@cmd /C "@echo ### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@cmd /C "@echo ### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files..."
	@if exist $(PRODUCT) $(RM) $(PRODUCT)
	$(RM) $(ALL_OBJS)
	$(ECHO) "### Deleted all derived files."


