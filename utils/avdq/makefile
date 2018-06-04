CC=g++

TARGET := avdq
source_dirs := source ../../source
sources := avdq.cpp math_addons.cpp atmosphere.cpp avdheatflux.cpp magnitude.cpp
includes := source ../../include
include_dirs := $(foreach d, $(includes), -I$d)
object_cpp_files := $(notdir $(sources) )
object_files := $(object_cpp_files:.cpp=.o)


VPATH := $(source_dirs)

all: $(object_files)
	$(CC) $^ -lm -o $(TARGET)
	rm *.o *.d
%.o: %.cpp
	$(CC) -c $(include_dirs) -Wno-deprecated $< -MD
	
include $(wildcard *.d)

clean:
	rm *.o *.d
