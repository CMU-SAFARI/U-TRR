program_NAME := TRRAnalyzer
program_CXX_SRCS := TRRAnalyzer.cpp $(wildcard ${DRAM_BENDER_ROOT}/sources/api/*.c) $(wildcard ${DRAM_BENDER_ROOT}/sources/api/*.cpp)
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_CXX_OBJS := ${program_CXX_OBJS:.c=.o}
program_OBJS := $(program_CXX_OBJS)
program_INCLUDE_DIRS := ${DRAM_BENDER_ROOT}/sources/api ../
program_LIBRARIES := pthread boost_program_options boost_filesystem boost_system
CPPFLAGS += -g -O3 -std=c++11

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

CC=g++

.PHONY: all clean distclean

all: $(program_OBJS)
	$(CC) $(CPPFLAGS) $(program_OBJS) -o $(program_NAME) $(LDFLAGS)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean
