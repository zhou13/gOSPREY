# Generated automatically from Makefile.in by configure.
# File         : Makefile
# Author       : Sung-Hoon Ko, Xinying Li
# Created      : Thu Apr 1998
# Revision     : ver 1.1
# Updated      : $Date: 2002/01/18 23:09:32 $
# Copyright: Northeast Parallel Architectures Center at Syracuse University 1998

SHELL       = /bin/sh


srcdir = .
top_srcdir = ../../..

VPATH = .:${srcdir}


JDK = /usr/java/latest

MPI_INC = -I/usr/pkg/mpich2-1.0.7-x86_64/include
MPI_LIB = -L/usr/pkg/mpich2-1.0.7-x86_64/lib

CCLIBKEY = 

INCLUDES = -I$(JDK)/include -I$(JDK)/include/linux -I$(MPI_INC)

CLASSES = ../../../lib/classes
CLASSPATH = .:$(CLASSES)

LIB = ../../../lib

CDIR = ../../../src/C

TITLE = mpi
PFX = $(TITLE)_
PKG = $(TITLE).
DIR = $(TITLE)/

RM = rm
CP = cp
MV = mv

JC = $(JDK)/bin/javac
JH = $(JDK)/bin/javah
JFLAG = -jni
CFLAGS = -c


.SUFFIXES:
.SUFFIXES: .java .h .class


JAVAFILES= \
     MPI.java           \
     CartParms.java     \
     Cartcomm.java      \
     Comm.java          \
     Datatype.java      \
     Errhandler.java    \
     GraphParms.java    \
     Graphcomm.java     \
     Group.java         \
     Intercomm.java     \
     Intracomm.java     \
     Op.java            \
     Prequest.java      \
     Request.java       \
     ShiftParms.java    \
     Status.java        \
     User_function.java 

JAVA_H= $(CDIR)/$(PFX)MPI.h           \
        $(CDIR)/$(PFX)CartParms.h     \
        $(CDIR)/$(PFX)Cartcomm.h      \
        $(CDIR)/$(PFX)Comm.h          \
        $(CDIR)/$(PFX)Datatype.h      \
        $(CDIR)/$(PFX)Errhandler.h    \
        $(CDIR)/$(PFX)GraphParms.h    \
        $(CDIR)/$(PFX)Graphcomm.h     \
        $(CDIR)/$(PFX)Group.h         \
        $(CDIR)/$(PFX)Intercomm.h     \
        $(CDIR)/$(PFX)Intracomm.h     \
        $(CDIR)/$(PFX)Op.h            \
        $(CDIR)/$(PFX)Prequest.h      \
        $(CDIR)/$(PFX)Request.h       \
        $(CDIR)/$(PFX)ShiftParms.h    \
        $(CDIR)/$(PFX)Status.h        \
        $(CDIR)/$(PFX)User_function.h 



all: MPIJAVACLASS $(JAVA_H)


MPIJAVACLASS : $(JAVAFILES)
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JC) -d $(CLASSES)   *.java 

$(CDIR)/$(PFX)MPI.h : $(CLASSES)/$(DIR)MPI.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)MPI
$(CDIR)/$(PFX)CartParms.h : $(CLASSES)/$(DIR)CartParms.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)CartParms
$(CDIR)/$(PFX)Cartcomm.h  : $(CLASSES)/$(DIR)Cartcomm.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Cartcomm 	
$(CDIR)/$(PFX)Comm.h : $(CLASSES)/$(DIR)Comm.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Comm 
$(CDIR)/$(PFX)Datatype.h : $(CLASSES)/$(DIR)Datatype.class 
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Datatype 
$(CDIR)/$(PFX)Errhandler.h : $(CLASSES)/$(DIR)Errhandler.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Errhandler
$(CDIR)/$(PFX)GraphParms.h : $(CLASSES)/$(DIR)GraphParms.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)GraphParms
$(CDIR)/$(PFX)Graphcomm.h : $(CLASSES)/$(DIR)Graphcomm.class  
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Graphcomm 
$(CDIR)/$(PFX)Group.h   : $(CLASSES)/$(DIR)Group.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Group
$(CDIR)/$(PFX)Intercomm.h : $(CLASSES)/$(DIR)Intercomm.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Intercomm
$(CDIR)/$(PFX)Intracomm.h : $(CLASSES)/$(DIR)Intracomm.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Intracomm
$(CDIR)/$(PFX)Op.h     : $(CLASSES)/$(DIR)Op.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Op
$(CDIR)/$(PFX)Prequest.h : $(CLASSES)/$(DIR)Prequest.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Prequest
$(CDIR)/$(PFX)Request.h   : $(CLASSES)/$(DIR)Request.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Request
$(CDIR)/$(PFX)ShiftParms.h : $(CLASSES)/$(DIR)ShiftParms.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)ShiftParms
$(CDIR)/$(PFX)Status.h    : $(CLASSES)/$(DIR)Status.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)Status
$(CDIR)/$(PFX)User_function.h : $(CLASSES)/$(DIR)User_function.class
	CLASSPATH=$(CLASSPATH) ; export CLASSPATH ; \
	$(JH) $(JFLAG) -d $(CDIR) $(PKG)User_function


clean :
	$(RM) -f $(CLASSES)/$(DIR)*.class

# Things to do
#
#   Find out why we are setting CLASSPATH to itself all the time...
#

