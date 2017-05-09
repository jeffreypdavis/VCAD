#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Facet.o \
	${OBJECTDIR}/Facet_2D.o \
	${OBJECTDIR}/Facet_3D.o \
	${OBJECTDIR}/Intersect_Meshes_2D.o \
	${OBJECTDIR}/Intersect_Meshes_3D.o \
	${OBJECTDIR}/Mesh_2D.o \
	${OBJECTDIR}/Mesh_3D.o \
	${OBJECTDIR}/Point_2D.o \
	${OBJECTDIR}/Point_3D.o \
	${OBJECTDIR}/Simplify_Mesh_2D.o \
	${OBJECTDIR}/Simplify_Mesh_3D.o \
	${OBJECTDIR}/VSCAD_Error.o \
	${OBJECTDIR}/Valid_Mesh_2D.o \
	${OBJECTDIR}/Valid_Mesh_3D.o \
	${OBJECTDIR}/Vector_2D.o \
	${OBJECTDIR}/Vector_3D.o \
	${OBJECTDIR}/shapes.o \
	${OBJECTDIR}/stl.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVCAD_lib.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVCAD_lib.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libVCAD_lib.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/Facet.o: Facet.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Facet.o Facet.cpp

${OBJECTDIR}/Facet_2D.o: Facet_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Facet_2D.o Facet_2D.cpp

${OBJECTDIR}/Facet_3D.o: Facet_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Facet_3D.o Facet_3D.cpp

${OBJECTDIR}/Intersect_Meshes_2D.o: Intersect_Meshes_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Intersect_Meshes_2D.o Intersect_Meshes_2D.cpp

${OBJECTDIR}/Intersect_Meshes_3D.o: Intersect_Meshes_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Intersect_Meshes_3D.o Intersect_Meshes_3D.cpp

${OBJECTDIR}/Mesh_2D.o: Mesh_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Mesh_2D.o Mesh_2D.cpp

${OBJECTDIR}/Mesh_3D.o: Mesh_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Mesh_3D.o Mesh_3D.cpp

${OBJECTDIR}/Point_2D.o: Point_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Point_2D.o Point_2D.cpp

${OBJECTDIR}/Point_3D.o: Point_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Point_3D.o Point_3D.cpp

${OBJECTDIR}/Simplify_Mesh_2D.o: Simplify_Mesh_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Simplify_Mesh_2D.o Simplify_Mesh_2D.cpp

${OBJECTDIR}/Simplify_Mesh_3D.o: Simplify_Mesh_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Simplify_Mesh_3D.o Simplify_Mesh_3D.cpp

${OBJECTDIR}/VSCAD_Error.o: VSCAD_Error.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/VSCAD_Error.o VSCAD_Error.cpp

${OBJECTDIR}/Valid_Mesh_2D.o: Valid_Mesh_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Valid_Mesh_2D.o Valid_Mesh_2D.cpp

${OBJECTDIR}/Valid_Mesh_3D.o: Valid_Mesh_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Valid_Mesh_3D.o Valid_Mesh_3D.cpp

${OBJECTDIR}/Vector_2D.o: Vector_2D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Vector_2D.o Vector_2D.cpp

${OBJECTDIR}/Vector_3D.o: Vector_3D.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Vector_3D.o Vector_3D.cpp

${OBJECTDIR}/shapes.o: shapes.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/shapes.o shapes.cpp

${OBJECTDIR}/stl.o: stl.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/stl.o stl.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
