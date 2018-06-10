# PATH
EXECPATH = ./bin/

# source files
SOURCES = ./src/*.cpp

# executable name
EXECNAME = egn2

# parameters for g++ compiler
GCC_NAME    = g++
GCC_FLAGS   = -std=c++11
GCC_LIBS    = -pthread

.SUFFIXES: .cpp

exec: $(OBJECTS)
	$(GCC_NAME) $(SOURCES) $(GCC_LIBS) $(GCC_FLAGS) -o $(EXECPATH)$(EXECNAME)
clean:
	rm $(EXECNAME)
# END
