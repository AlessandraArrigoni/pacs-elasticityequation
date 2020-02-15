# Compiler
CXX = g++ -std=c++14

# Optimize flags
OPTFLAGS = -g
GETFEM_PATH = /home/alearrigoni/Documenti/getfem-5.3
INCLUDE = -I$(GETFEM_PATH)/include -I$(GETFEM_PATH)/src -I$(GETFEM_PATH)/src/gmm -I../include -I/usr/include
# Flags
CXXFLAGS = $(INCLUDE)  $(OPTFLAGS)

# Executable source
EXESRCS = *.cc

# Executable object file
EXEOBJS = $(EXESRCS:.cc = .o)

# Executable name
EXEC = LaplacianBasic

# Sources folder
FOLDER = ../src/

# Laptop
LIB_PATH =  /usr/lib

# Laptop
GETFEM_LIB =  /usr/local/lib/libgetfem.a

# Laptop
LDLIBS += $(GETFEM_LIB) -rdynamic /usr/lib/x86_64-linux-gnu/libqhull.so.7 /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/libblas.so.3


DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS

# Sources
SRCS = $(wildcard $(FOLDER)*.cc)

# Objects
OBJS = $(SRCS:.cc=.o)

# Headers
HEADERS = $(SRCS:.cc=.h)

# Name file of dependences
DEPEND = make.dep

.PHONY: all clean optimised debug test

all : $(DEPEND) $(OBJS) $(EXEOBJS)
	$(CXX) $(OPTFLAGS) -o $(EXEC) $(EXEOBJS) $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)

$(DEPEND) : $(SRCS) $(EXESRCS)
	$(CXX) -MM $(SRCS) $(EXESRCS) -MF $(DEPEND)  $(INCLUDE)

-include $(DEPEND)

debug:
	@$(MAKE) --no-print-directory OPTFLAGS+="-DDEBUG" all

optimised :
	@$(MAKE) --no-print-directory OPTFLAGS="-O3 -DNDEBUG" all

test :
	@$(MAKE) --no-print-directory OPTFLAGS="-O3 -DTEST -DNDEBUG" all

clean :
	$(RM) $(EXEC) $(OBJS) $(EXECOBJS) outputData/*.mm output_vtk/*.vtk outputData/laplacianBasic_errors*
