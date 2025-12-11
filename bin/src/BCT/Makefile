CXX = g++
CXXFLAGS = -O3 -std=c++11
TARGET = ../../ctw-calc

all: $(TARGET)

$(TARGET): main_ctw.cpp Functions.h node.h
	$(CXX) $(CXXFLAGS) -o $(TARGET) main_ctw.cpp

clean:
	rm -f $(TARGET)
