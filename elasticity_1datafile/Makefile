# Compiler
CXX = g++

# Optimize flags
OPTFLAGS = -g #-O3
GETFEM_PATH = /home/sara/Scrivania/getfem-5.3
INCLUDE = -I$(GETFEM_PATH)/include -I$(GETFEM_PATH)/src -I$(GETFEM_PATH)/src/gmm -I./include -I/usr/include
# Flags
CXXFLAGS = $(INCLUDE)  $(OPTFLAGS)

# Executable source
EXESRCS = *.cc

# Executable object file
EXEOBJS = $(EXESRCS:.cc = .o)

# Executable name
EXEC = linElasticityIFace

# Sources folder
FOLDER = src/


# Laptop
LIB_PATH =  /usr/lib


# Laptop
LDLIBS =  /usr/local/lib/libgetfem.a

# Laptop
LDLIBS += $(GETFEM_LIB) -rdynamic    /usr/lib/x86_64-linux-gnu/libqhull.so.7 /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/libblas.so.3


DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS

# Sources
SRCS = $(wildcard $(FOLDER)*.cc)

# Objects
OBJS = $(SRCS:.cc=.o)

# Headers
HEADERS = $(SRCS:.cc=.h)

# Name file of dependences
DEPEND = make.dep

.PHONY: all clean

all : $(DEPEND) $(OBJS) $(EXEOBJS)
	$(CXX) $(OPTFLAGS) -o $(EXEC) $(EXEOBJS) $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)

decoupled  :  $(DEPEND) $(OBJS) main_decoupled.o
	$(CXX) $(OPTFLAGS) -o $(EXEC) main_decoupled.o $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)

coupled  :  $(DEPEND) $(OBJS) main_coupled.o
	$(CXX) $(OPTFLAGS) -o $(EXEC) main_coupled.o $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)

iterative  :  $(DEPEND) $(OBJS) main_iterative.o
	$(CXX) $(OPTFLAGS) -o $(EXEC) main_iterative.o $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)

test  :  $(DEPEND) $(OBJS) main_test.o
	$(CXX) $(OPTFLAGS) -o $(EXEC) main_test.o $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE)


$(DEPEND) : $(SRCS) $(EXESRCS)
	$(CXX) -MM $(SRCS) $(EXESRCS) -MF $(DEPEND)  $(INCLUDE)

-include $(DEPEND)

clean :
	-rm $(EXEC) $(OBJS) $(EXECOBJS)
