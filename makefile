# Compiler settings - Can be customized.
CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -fcompare-debug-second

# Google Test settings
GTEST_DIR = C:\Users\James\OneDrive\VSCode\ENGR570\Final project\-MaxwellSolver_NERS570\googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# Targets
all: main test

# Main application
main: src/FDTD.o src/Point.o src/CurrentLine.o src/Mesh.o src/main.o src/Source.o src/globals.o
	$(CXX) $(CXXFLAGS) -o main src/FDTD.o src/main.o src/Mesh.o src/Source.o src/Point.o src/CurrentLine.o src/globals.o

# Test
test: src/PointTest.o src/CurrentLineTest.o
	$(CXX) $(CXXFLAGS) -pthread -o test src/PointTest.o src/CurrentLineTest.o libgtest.a

# Object files
src/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp -o src/main.o

src/FDTD.o: src/FDTD.cpp include/FDTD.hpp
	$(CXX) $(CXXFLAGS) -c src/FDTD.cpp -o src/FDTD.o

src/Point.o: src/Point.cpp include/Point.hpp
	$(CXX) $(CXXFLAGS) -c src/Point.cpp -o src/Point.o

src/CurrentLine.o: src/CurrentLine.cpp include/CurrentLine.hpp
	$(CXX) $(CXXFLAGS) -c src/CurrentLine.cpp -o src/CurrentLine.o

src/Mesh.o: src/Mesh.cpp include/Mesh.hpp
	$(CXX) $(CXXFLAGS) -c src/Mesh.cpp -o src/Mesh.o

src/Source.o: src/Source.cpp include/Source.hpp
	$(CXX) $(CXXFLAGS) -c src/Source.cpp -o src/Source.o

src/globals.o: src/globals.cpp include/globals.hpp
	$(CXX) $(CXXFLAGS) -c src/globals.cpp -o src/globals.o

src/PointTest.o: tests/PointTest.cpp
	$(CXX) $(CXXFLAGS) -I$(GTEST_DIR)/include -c tests/PointTest.cpp -o src/PointTest.o

src/CurrentLineTest.o: tests/CurrentLineTest.cpp
	$(CXX) $(CXXFLAGS) -I$(GTEST_DIR)/include -c tests/CurrentLineTest.cpp -o src/CurrentLineTest.o

# Google Test
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CXXFLAGS) -I$(GTEST_DIR) -I$(GTEST_DIR)/include -pthread -c \
            $(GTEST_DIR)/src/gtest-all.cc

libgtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

# Clean
clean:
	rm -f src/*.o main test libgtest.a
