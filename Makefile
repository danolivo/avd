CC=g++
# CC=C:\GCC\BIN\g++
# CC=C:\MinGW\BIN\mingw32-g++

source_dirs := . source source/1DThermal
includes := include
include_dirs := $(foreach d, $(includes), -I$d)

wild := $(addsuffix /*.cpp, $(source_dirs) )
source_cpp_files := $(wildcard $(addsuffix /*.cpp, $(source_dirs) ) )
source_c_files := $(wildcard $(addsuffix /*.c, $(source_dirs) ) )
object_c_files := $(notdir $(source_c_files) )
object_cpp_files := $(notdir $(source_cpp_files) )
object_files := $(object_cpp_files:.cpp=.o) $(object_c_files:.c=.o)


VPATH := $(source_dirs)

all: $(object_files)
	$(CC) $^ -lm -o avd
	rm *.o *.d
%.o: %.cpp
	$(CC) -c $(include_dirs) -Wno-deprecated $< -MD
%.o: %.c
	$(CC) -c $(include_dirs) -Wno-deprecated $< -MD
	
include $(wildcard *.d)

clean:
	rm *.o *.d

