CC=g++

target := xeff_base
source_dirs := source ../../source
sources := source/main.cpp source/xeff_table.cpp ../source/math_addons.cpp
includes := source ../../include
include_dirs := $(foreach d, $(includes), -I$d)
object_cpp_files := $(notdir $(sources) )
object_files := $(object_cpp_files:.cpp=.o)


VPATH := $(source_dirs)

all: $(object_files)
	$(CC) $^ -lm -o $(target)
	rm *.o *.d
%.o: %.cpp
	$(CC) -c $(include_dirs) -Wno-deprecated $< -MD
	
include $(wildcard *.d)

clean:
	rm *.o *.d
