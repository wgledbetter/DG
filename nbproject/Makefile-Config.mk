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
CND_CONF=Config
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/pch.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-march=native -mtune=native -std=c++17 -fPIC -pthread
CXXFLAGS=-march=native -mtune=native -std=c++17 -fPIC -pthread

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/home/wgledbetter/Install/intel/mkl/lib/intel64 -L../../Programs/SuiteSparse/SuiteSparse/lib -Wl,-rpath,'../../Install/intel/mkl/lib/intel64' -Wl,-rpath,'../../Programs/SuiteSparse/SuiteSparse/lib' -lmkl_core -lmkl_sequential -lmkl_intel_lp64 -lspqr -lamd -lcholmod

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dg

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dg: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dg ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DEIGEN_USE_MKL_ALL -I../../Install/intel/mkl/include -I../../Programs/SuiteSparse/SuiteSparse/include -I../../Install/Eigen/eigen-git-mirror -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/pch.o: src/pch.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DEIGEN_USE_MKL_ALL -I../../Install/intel/mkl/include -I../../Programs/SuiteSparse/SuiteSparse/include -I../../Install/Eigen/eigen-git-mirror -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pch.o src/pch.cpp

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
